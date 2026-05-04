/// Coupling modules: cross-species interactions.
///
/// Each coupling knows how to apply an interaction between species.
/// The engine composes couplings with free evolution in a Strang
/// splitting loop.
pub mod poisson;
