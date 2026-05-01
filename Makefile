.PHONY: lint test build clean reset publish

# Code quality
lint:
	cargo fmt
	cargo clippy --fix --allow-dirty

test:
	cargo test

# Build and release
build:
	cargo build --release --features vis
	@echo ""
	@echo "Run with: ./target/release/hermes"

publish:
	@VERSION=$$(grep '^version' Cargo.toml | head -1 | sed 's/.*"\(.*\)"/\1/'); \
	if git rev-parse "v$$VERSION" >/dev/null 2>&1; then \
		echo "Error: Tag v$$VERSION already exists"; \
		echo "Update the version in Cargo.toml first"; \
		exit 1; \
	fi; \
	echo "Publishing version $$VERSION"; \
	git tag "v$$VERSION" && \
	git push origin main --tags

# Cleanup
clean:
	cargo clean
	find . -type f -name ".DS_Store" -delete 2>/dev/null || true

reset: clean build
