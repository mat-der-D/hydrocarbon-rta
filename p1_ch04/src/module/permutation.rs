use itertools::Itertools;
use std::collections::HashSet;

use super::matrix::AdjacencyMatrix;

pub fn calc_orbit<const N: usize>(mat: AdjacencyMatrix<N>) -> HashSet<AdjacencyMatrix<N>> {
    let mut orbit = HashSet::new();
    for perm_vec in (0..N).permutations(N) {
        let mut perm = [0; N];
        perm.copy_from_slice(&perm_vec);
        orbit.insert(mat.permute_indices(&perm));
    }
    orbit
}
