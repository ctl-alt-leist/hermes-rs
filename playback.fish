#!/usr/bin/env fish

# Play back saved hermes snapshots from a directory.
#
# Usage:
#   ./playback.fish cosmic-web-pm
#   ./playback.fish scenes/cosmic-web-pm
#   ./playback.fish scenes/cosmic-web-pm --fps 60
#   ./playback.fish scenes/cosmic-web-pm --record cosmic-web.gif

if test (count $argv) -lt 1
    echo "Usage: playback.fish <directory> [hermes flags...]"
    echo ""
    echo "Examples:"
    echo "  playback.fish cosmic-web-pm"
    echo "  playback.fish scenes/cosmic-web-pm --fps 60"
    echo "  playback.fish scenes/cosmic-web-pm --record output.gif"
    exit 1
end

set dir $argv[1]
set rest $argv[2..]

# Resolve bare names to scenes/<name>
if not string match -q '*/*' $dir
    set dir "scenes/$dir"
end

if not test -d $dir
    echo "Directory not found: $dir"
    exit 1
end

cargo run --release --features vis -- --playback $dir $rest
