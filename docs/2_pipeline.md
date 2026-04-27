# Pipeline Architecture

The simulation, disk I/O, and visualization run on independent threads connected by bounded channels. The simulation never blocks on consumers; consumers never block the simulation. Adding a new consumer (remote viewer, metrics exporter) is: write one function that reads from a channel, connect it to the router.

## Topology

```
Sim Thread ──ch(8)──> Router ──ch(16)──> Disk Writer
                         └────ch(4)──> Precompute ──ch(4)──> Main Thread (Viewer)
```

Each arrow is a `mpsc::sync_channel` with the capacity shown. The simulation produces `Arc<Snapshot>` into its output channel. The router fans out to consumers by cloning the `Arc` (only the reference count is incremented -- no particle data is copied). Consumers that fall behind have frames silently dropped via `try_send`.

## Threads by Mode

| Mode | Sim | Router | Disk | Precompute | Main |
|------|-----|--------|------|------------|------|
| `hermes --save` | spawned | spawned | spawned | -- | blocks on join |
| `hermes --live` | spawned | spawned | optional | spawned | viewer loop |
| `hermes --playback` | -- | -- | -- | loader | viewer loop |
| `hermes --record` | -- | -- | -- | -- | sequential encode |

The simulation always runs on a spawned thread. The main thread always owns the event loop (macOS requires the window event loop on the main thread). In headless mode, the main thread simply blocks on `sim_handle.join()`.

## Key Types

### PipelineMessage

```rust
pub enum PipelineMessage {
    Snapshot(Arc<Snapshot>),
    Done,
}
```

The `Arc<Snapshot>` wrapping enables zero-copy fan-out. A `Snapshot` for 32K particles is ~500 KB; cloning it for each consumer would double memory traffic. With `Arc`, the router just increments a reference count.

### DisplayFrame

```rust
pub struct DisplayFrame {
    pub positions: Vec<[f32; 3]>,
    pub colors: Vec<[f32; 3]>,
    pub step: usize,
    pub scale_factor: f64,
}
```

The precompute thread converts morphis `Vector<3>` positions and velocity-based colors into flat f32 arrays. The viewer thread never touches morphis types -- it just draws from flat arrays.

### SnapshotSender

```rust
pub struct SnapshotSender {
    tx: SyncSender<PipelineMessage>,
}
```

The simulation's output interface. `send()` uses `try_send` (non-blocking, drops on full). `done()` uses blocking `send` to ensure the shutdown signal reaches the router.

## Shutdown Ordering

1. Simulation finishes its step loop, calls `sender.done()`
2. Router receives `Done`, forwards to all consumers via blocking send, exits
3. Disk writer receives `Done`, prints summary, exits
4. Precompute thread receives `Done`, sends `ViewerMessage::Done`, exits
5. Main thread viewer sees `Done` or user closes window, returns
6. `main()` joins all handles

If the user closes the viewer window early, the simulation and disk writer continue to completion. This is correct: closing the window should not abort a running simulation or discard unsaved data.

## Playback

The playback viewer uses a loader thread that reads snapshots from disk, precomputes `DisplayFrame`s, and sends them via a bounded channel. The main thread collects all frames (with a progress bar) before starting the render loop. This ensures smooth playback with zero I/O during rendering.

## Performance

Measured at 32K particles, 50 steps:

| Path | Time | Notes |
|------|------|-------|
| Observer + FileObserver (old) | 17.5 s | Sync writes block simulation |
| Pipeline + disk writer (new) | 3.2 s | Writes on separate thread |

The pipeline is **5.5x faster** with disk saving enabled because the simulation never waits for I/O. Without saving, the overhead of `Arc::new(snapshot)` is negligible compared to the simulation step cost.
