/// Error types for hermes.
#[derive(Debug, thiserror::Error)]
pub enum HermesError {
    #[error("configuration error: {0}")]
    Config(String),

    #[error("cosmology error: {0}")]
    Cosmology(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),
}
