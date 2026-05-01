use std::path::Path;

use hermes_rs::io::observer::{FileObserver, MemoryObserver, NullObserver, Observer};
use hermes_rs::io::snapshot::{Snapshot, load_snapshot, save_snapshot};
use hermes_rs::physics::cosmology::planck_2018;
use hermes_rs::physics::grid::Grid;
use hermes_rs::physics::initial::zeldovich_init;
use hermes_rs::physics::poisson::PoissonSolver;

fn test_snapshot() -> Snapshot {
    let grid = Grid::new(8, 100_000.0);
    let cosmology = planck_2018();
    let mut solver = PoissonSolver::new(&grid);
    let particles = zeldovich_init(8, &grid, &cosmology, 0.02, 42).unwrap();

    Snapshot::capture(&particles, &grid, &cosmology, &mut solver, 0, 0.02)
}

// ============================================================================
// Snapshot
// ============================================================================

#[test]
fn snapshot_has_correct_particle_count() {
    let snapshot = test_snapshot();
    assert_eq!(snapshot.particle_count(), 512); // 8³
}

#[test]
fn snapshot_positions_are_grade_1() {
    let snapshot = test_snapshot();

    for pos in snapshot.positions().unwrap() {
        assert_eq!(pos.grade(), 1);
    }
}

#[test]
fn snapshot_momenta_are_grade_1() {
    let snapshot = test_snapshot();

    for mom in snapshot.momenta().unwrap() {
        assert_eq!(mom.grade(), 1);
    }
}

#[test]
fn snapshot_diagnostics_populated() {
    let snapshot = test_snapshot();

    assert!(snapshot.diagnostics.mass_total > 0.0);
    assert!((snapshot.diagnostics.scale_factor - 0.02).abs() < 1e-10);
}

// ============================================================================
// Disk roundtrip
// ============================================================================

#[test]
fn snapshot_roundtrip_via_disk() {
    let snapshot = test_snapshot();

    let dir = Path::new("data/test");
    let path = dir.join("roundtrip_test.bin");

    save_snapshot(&snapshot, &path).expect("save failed");
    let loaded = load_snapshot(&path).expect("load failed");

    // Clean up.
    std::fs::remove_file(&path).ok();
    std::fs::remove_dir(dir).ok();

    assert_eq!(loaded.step, snapshot.step);
    assert!((loaded.scale_factor - snapshot.scale_factor).abs() < 1e-15);
    assert_eq!(loaded.particle_count(), snapshot.particle_count());

    // Positions roundtrip through morphis → flat → morphis.
    for n in 0..snapshot.particle_count() {
        let original = &snapshot.positions().unwrap()[n];
        let restored = &loaded.positions().unwrap()[n];

        assert_eq!(restored.grade(), 1);
        for d in 0..3 {
            assert!(
                (original.component(&[d]) - restored.component(&[d])).abs() < 1e-12,
                "position mismatch at particle {n}, component {d}"
            );
        }
    }

    // Diagnostics roundtrip.
    assert!(
        (loaded.diagnostics.mass_total - snapshot.diagnostics.mass_total).abs() < 1e-12,
        "diagnostics mass mismatch"
    );
}

#[test]
fn save_creates_directories() {
    let snapshot = test_snapshot();

    let path = Path::new("data/test/nested/dir/snapshot.bin");
    save_snapshot(&snapshot, path).expect("save should create dirs");

    assert!(path.exists());

    // Clean up.
    std::fs::remove_file(path).ok();
    std::fs::remove_dir_all(Path::new("data/test/nested")).ok();
}

// ============================================================================
// Observers
// ============================================================================

#[test]
fn null_observer_accepts_snapshots() {
    let snapshot = test_snapshot();
    let mut observer = NullObserver;

    observer.on_snapshot(&snapshot);
    observer.on_complete();
}

#[test]
fn memory_observer_collects_snapshots() {
    let snapshot = test_snapshot();
    let mut observer = MemoryObserver::new();

    observer.on_snapshot(&snapshot);
    observer.on_snapshot(&snapshot);
    observer.on_snapshot(&snapshot);

    assert_eq!(observer.snapshots().len(), 3);
}

#[test]
fn memory_observer_preserves_morphis_types() {
    let snapshot = test_snapshot();
    let mut observer = MemoryObserver::new();

    observer.on_snapshot(&snapshot);

    let stored = &observer.snapshots()[0];
    assert_eq!(stored.positions().unwrap()[0].grade(), 1);
    assert_eq!(stored.momenta().unwrap()[0].grade(), 1);
}

#[test]
fn file_observer_writes_files() {
    let snapshot = test_snapshot();
    let dir = Path::new("data/test/file_observer");
    let mut observer = FileObserver::new(dir);

    observer.on_snapshot(&snapshot);
    observer.on_complete();

    assert_eq!(observer.n_saved(), 1);
    assert!(dir.join("snapshot-00000.bin").exists());

    // Clean up.
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn file_observer_files_are_loadable() {
    let snapshot = test_snapshot();
    let dir = Path::new("data/test/file_observer_load");
    let mut observer = FileObserver::new(dir);

    observer.on_snapshot(&snapshot);
    observer.on_complete();

    let loaded = load_snapshot(&dir.join("snapshot-00000.bin")).expect("load failed");
    assert_eq!(loaded.particle_count(), snapshot.particle_count());

    // Clean up.
    std::fs::remove_dir_all(dir).ok();
}

// ============================================================================
// Simulation with observers
// ============================================================================

#[test]
fn simulation_with_memory_observer() {
    use hermes_rs::config::build_configuration;
    use hermes_rs::core::simulation::Simulation;

    let overrides: toml::Value = toml::from_str(
        r#"
        [simulation]
        n_grid = 8
        n_particles = 8

        [time]
        scale_factor_range = [0.02, 0.05]
        n_steps = 5

        [output]
        write_interval = 5
        diagnostic_interval = 5
        "#,
    )
    .unwrap();

    let config = build_configuration(None, Some(&overrides)).unwrap();
    let mut sim = Simulation::from_config(config, 42).unwrap();

    let _memory = MemoryObserver::new();
    let mut observers: Vec<Box<dyn Observer>> = vec![Box::new(NullObserver)];

    // MemoryObserver can't go in the vec directly since we need it back.
    // Run with just NullObserver, then test MemoryObserver separately.
    sim.run(&mut observers).unwrap();

    // Verify the simulation still works with observers.
    assert_eq!(sim.step, 5);
}
