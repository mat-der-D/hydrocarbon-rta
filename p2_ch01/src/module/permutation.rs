use itertools::Itertools;
use rustc_hash::FxHashSet;

use super::matrix::AdjacencyMatrix;

pub fn calc_orbit<const N: usize>(mat: AdjacencyMatrix<N>) -> FxHashSet<AdjacencyMatrix<N>> {
    let mut orbit = FxHashSet::default();
    for perm_vec in (0..N).permutations(N) {
        let mut perm = [0; N];
        perm.copy_from_slice(&perm_vec);
        orbit.insert(mat.permute_indices(&perm));
    }
    orbit
}
