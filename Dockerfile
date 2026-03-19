# Build stage - use musl target for static binary
FROM docker.io/rust:1.94-alpine AS builder
WORKDIR /build

# Copy dependency manifests and source code
COPY Cargo.toml Cargo.lock* ./
COPY src ./src

# Build static release binary
RUN cargo build --release

# Runtime stage - use samtools as base
FROM docker.io/staphb/samtools:1.23

# Copy kineticsquash binary from builder
COPY --from=builder /build/target/release/kineticsquash /usr/local/bin/kineticsquash

ENTRYPOINT ["/usr/local/bin/kineticsquash"]
