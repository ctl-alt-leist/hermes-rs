#!/usr/bin/env fish

# Run a hermes simulation from a scene TOML.
#
# Usage:
#   ./run.fish cosmic-web-pm
#   ./run.fish scenes/cosmic-web-ft.toml
#   ./run.fish cosmic-web-pm --live --save
#   ./run.fish cosmic-web-pm --steps 100 --seed 7

if test (count $argv) -lt 1
    echo "Usage: run.fish <scene> [hermes flags...]"
    echo ""
    echo "Examples:"
    echo "  run.fish cosmic-web-pm"
    echo "  run.fish cosmic-web-pm --live --save"
    echo "  run.fish scenes/galaxy-group-ft.toml --steps 200"
    exit 1
end

set scene $argv[1]
set rest $argv[2..]

# Resolve bare names to scenes/<name>.toml
if not string match -q '*/*' $scene
    set scene "scenes/$scene"
end
if not string match -q '*.toml' $scene
    set scene "$scene.toml"
end

if not test -f $scene
    echo "Scene not found: $scene"
    exit 1
end

cargo run --release --features vis -- --scene $scene $rest
