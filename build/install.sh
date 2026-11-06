#!/bin/sh
set -eu
app_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
command_dir="$HOME/.local/bin"
mkdir -p "$command_dir"
for command in myporo mycfs; do
    test -x "$app_dir/bin/$command"
    if [ -e "$command_dir/$command" ] && [ ! -L "$command_dir/$command" ]; then
        echo "A file already exists: $command_dir/$command" >&2
        exit 1
    fi
    ln -sfn "$app_dir/bin/$command" "$command_dir/$command"
done
printf 'Commands installed in %s\n' "$command_dir"
case ":$PATH:" in
    *":$command_dir:"*) ;;
    *) printf 'Add this line to your shell configuration, then open a new terminal:\nexport PATH="$HOME/.local/bin:$PATH"\n' ;;
esac
