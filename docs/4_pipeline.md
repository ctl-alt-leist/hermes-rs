# Pipeline Architecture

The simulation, disk I/O, and visualization run on independent threads connected by bounded channels. The simulation never blocks on droppable consumers (viewer); non-droppable consumers (disk writer) apply backpressure if they fall behind. Adding a new consumer is: write one function that reads from a channel, connect it to the router.

## Topology

```
Sim Thread ──ch(512)──> Router ──ch(512)──> Disk Writer (non-droppable)
                           └────ch(4)──> Precompute ──ch(4)──> Main Thread (Viewer)
```

Each arrow is a `mpsc::sync_channel` with the capacity shown. The simulation produces `Arc<Snapshot>` into its output channel. The router fans out to consumers by cloning the `Arc` (only the reference count is incremented — no particle data is copied).

Consumer channels are configured as droppable or non-droppable via `ConsumerConfig`. The disk writer is non-droppable (blocking send — every snapshot reaches disk). The viewer precompute is droppable (`try_send` — frames are silently dropped when the channel is full, keeping the simulation running at full speed).

## Snapshot Gating

Not every simulation step produces a snapshot for the pipeline. The `write_interval` configuration controls which steps are sent to the disk writer. The simulation only captures and sends a snapshot when the step is a multiple of `write_interval` or is the final step. This reduces I/O without affecting the physics.

For resumed simulations, the step numbering continues from the last snapshot's step number, and the initial duplicate frame is skipped.

## Threads by Mode

| Mode | Sim | Router | Disk | Precompute | Main |
|------|-----|--------|------|------------|------|
| `cargo run -- --save` | spawned | spawned | spawned | -- | blocks on join |
| `cargo run -- --live` | spawned | spawned | optional | spawned | viewer loop |
| `cargo run -- --playback` | -- | -- | -- | loader | viewer loop |
| `cargo run -- --record` | -- | -- | -- | -- | sequential encode |
| `cargo run -- --resume` | spawned | spawned | spawned | -- | blocks on join |

The simulation always runs on a spawned thread. The main thread always owns the event loop (macOS requires the window event loop on the main thread). In headless and resume modes, the main thread blocks on `sim_handle.join()`.

## Rendering

The viewer uses kiss3d's `State` trait to support both particle and field rendering in the same window. The `ViewerState` dispatches on `RenderMode`:

- **Points**: Particle snapshots are rendered via `window.draw_point()` through the built-in kiss3d point renderer.
- **Volumetric**: Field snapshots are rendered via a custom `VolumetricRenderer` that draws point sprites with additive blending and Gaussian falloff. Each grid cell is a soft, semi-transparent blob; overlapping blobs accumulate brightness to produce a smooth volumetric appearance.

The `VolumetricRenderer` is returned from `cameras_and_effect_and_renderer()` on every frame but only does work when volumetric points have been queued. For particle frames it returns immediately with no GL state changes.

### Playback Controls

The playback viewer supports keyboard interaction:

| Key | Action |
|-----|--------|
| Space | Play / pause |
| Left / Right | Single frame step (auto-pauses) |
| Up / Down | Jump 10% forward / back (auto-pauses) |
| Home / End | Jump to first / last frame |

## Key Types

### PipelineMessage

```rust
pub enum PipelineMessage {
    Snapshot(Arc<Snapshot>),
    Done,
}
```

The `Arc<Snapshot>` wrapping enables zero-copy fan-out.

### DisplayFrame

```rust
pub struct DisplayFrame {
    pub positions: Vec<[f32; 3]>,
    pub colors: Vec<[f32; 3]>,
    pub render_mode: RenderMode,
    pub step: usize,
    pub scale_factor: f64,
}
```

The precompute thread converts morphis `Vector<3>` positions and velocity-based colors into flat f32 arrays, and tags the frame with its render mode. The viewer thread never touches morphis types.

### ConsumerConfig

```rust
pub struct ConsumerConfig {
    pub tx: SyncSender<PipelineMessage>,
    pub droppable: bool,
}
```

Controls whether the router blocks (disk) or drops (viewer) when the consumer's channel is full.

## Shutdown Ordering

1. Simulation finishes its step loop, calls `sender.done()`
2. Router receives `Done`, forwards to all consumers via blocking send, exits
3. Disk writer receives `Done`, exits
4. Precompute thread receives `Done`, sends `ViewerMessage::Done`, exits
5. Main thread viewer sees `Done` or user closes window, returns
6. `main()` joins all handles

Closing the viewer window early does not abort the simulation or discard unsaved data.

## Progress Bars

Headless and resume modes display two parallel progress bars (via indicatif `MultiProgress`): the simulation progress (cyan, with redshift/scale factor) and the disk writer progress (green, with snapshot count). The bars are vertically aligned.
