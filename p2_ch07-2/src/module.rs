mod dehydrogenation;
mod hot_bit_iter;
mod matrix;
mod parallel;
mod permutation;

pub use matrix::AdjacencyTwoBitsMatrix;
pub use parallel::{create_feat2skeletons, dehydrogenate_feat2skeletons};
pub use permutation::PermutationStore;
