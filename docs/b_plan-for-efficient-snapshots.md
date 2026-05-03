
CLAUDE: We'll spend more time focusing on the efficiency aspect of this only after we've unified the engine framework. However, as we go, feel free to recognize efficiencies we could take, or implement anything here that makes sense. Otherwise don't focus too much on this.

# Plan for Efficient Snapshots

## The Problem

At 64^3 particles, each snapshot is ~1.5 MB (bincode-serialized positions + momenta + metadata). Saving every step of a 300-step run produces 300 files totaling ~450 MB. The simulation itself takes ~2 minutes; the disk I/O adds significant time on top of that, and at higher resolutions the I/O will dominate.

The current architecture: one disk writer thread receives snapshots via a bounded channel and writes them sequentially. The simulation blocks if the channel fills (non-droppable consumer). This is correct but slow — the simulation can stall waiting for disk.

## Where the Bottleneck Might Be

Three candidates, and we don't yet know which dominates:

**1. Serialization (CPU-bound).** Bincode walks the `Snapshot` struct, copies morphis `Vector<3>` objects into flat `[f64; 3]` arrays, and serializes them. For 262,144 particles, this involves ~1.5M f64 copies plus bincode framing. This happens on the disk writer thread, blocking the next snapshot from being processed.

**2. Disk write (I/O-bound).** Writing 1.5 MB per file to SSD. Modern NVMe SSDs sustain 1-3 GB/s sequential write, so a single 1.5 MB write should take <1 ms. But file creation overhead (metadata, fsync, directory update) adds latency per file. With 300 files, the per-file overhead may matter more than the raw bytes.

**3. Channel backpressure.** The bounded channel (capacity 512) between router and disk writer. If serialization + write takes longer than the simulation step, the channel fills and the simulation stalls. The simulation step for 64^3 PM takes ~550 ms; if disk write + serialize takes >550 ms per snapshot, the simulation blocks.

## What to Measure

Before choosing a solution, benchmark the actual bottleneck on the target machine:

1. **Serialization time alone.** Serialize a snapshot to a `Vec<u8>` in memory (no disk write). Measure wall time. This isolates the CPU cost.

2. **Disk write time alone.** Write a pre-serialized `Vec<u8>` of the right size to a file. Measure wall time including file creation. Repeat for N files to capture per-file overhead.

3. **End-to-end pipeline throughput.** Run the simulation with the disk writer and measure how much the wall time exceeds the no-save case. The difference is the I/O overhead.

4. **Channel utilization.** Instrument the channel to count how often `send` blocks. If it never blocks, the I/O is keeping up and optimization isn't needed. If it blocks frequently, the writer is the bottleneck.

## Possible Solutions

### A. Reduce snapshot count

The simplest fix. `write_interval = 10` saves every 10th step, reducing I/O by 10x. The current default is 1 (every step). For playback at 30 fps, 30 snapshots cover 1 second of animation — 300 snapshots is 10 seconds of smooth playback. Writing every 3rd step would still give smooth playback and cut I/O by 3x.

Trade-off: temporal resolution. For analysis, every step matters. For visualization, every 3rd-5th step is fine.

### B. Parallelize serialization

If serialization is the bottleneck, move it off the disk writer thread. The simulation produces `Arc<Snapshot>` (cheap reference count). A pool of N serializer threads each take a snapshot, serialize it to `Vec<u8>`, and pass the bytes to the single disk writer. The writer just does `fs::write` — fast.

This separates CPU work (serialization, parallelizable) from I/O work (writing, sequential on one thread).

### C. Parallelize disk writes

If per-file overhead is the bottleneck, multiple writer threads each write to different files. SSDs handle parallel writes well (no seek penalty). Use a multi-consumer channel (crossbeam) or `Arc<Mutex<Receiver>>`.

Caution: too many parallel writes can cause write amplification on SSDs. 2-4 writer threads is probably the sweet spot.

### D. Batch writes

Instead of one file per snapshot, batch multiple snapshots into a single large file. Write a stream of serialized snapshots with a simple index (offset table at the start or end). This amortizes file creation overhead across many snapshots.

Trade-off: individual snapshots aren't independently loadable without reading the index. Playback and resume need to parse the batch format.

### E. Memory-mapped I/O

Use `mmap` to map the output file into virtual memory and write snapshot data directly. The OS handles the actual disk writes asynchronously. This can improve throughput for large sequential writes.

### F. Compression

Compress each snapshot (lz4, zstd). Particle data has spatial coherence — nearby particles have similar positions — so compression ratios of 2-5x are plausible. This trades CPU for disk bandwidth, which is favorable if disk is the bottleneck.

### G. Binary format optimization

The current format round-trips through morphis `Vector<3>` → `components_from_vector` → `[f64; 3]` → bincode. If we store the raw `Array2<f64>` (shape [3, N]) directly as a contiguous byte slice, serialization becomes a single memcpy. This eliminates the per-particle conversion loop.

## Recommended Investigation Order

1. Run the benchmark (A/B/C above under "What to Measure") to identify the actual bottleneck
2. If serialization dominates: try G (raw byte format) first, then B (parallel serialize)
3. If disk dominates: try A (write_interval) first, then D (batch writes)
4. If neither dominates (channel never blocks): the I/O is already fast enough, optimize elsewhere

## Relationship to Snapshot Format

The unified engine plan (see `unified-engine.md`) calls for a self-describing snapshot format with metadata headers. Any format change for efficiency should be designed alongside that — one migration, not two. A header + raw binary body is both self-describing and fast to serialize.
