use std::hash::Hash;

use itertools::Itertools;
use rustc_hash::FxHashSet;

pub trait Permutable<const N: usize> {
    fn permute_by(&self, perm: &[usize; N]) -> Self;
}

pub fn calc_orbit<const N: usize, P>(pable: P) -> FxHashSet<P>
where
    P: Permutable<N> + Eq + Hash,
{
    let mut orbit = FxHashSet::default();
    for perm_vec in (0..N).permutations(N) {
        let mut perm = [0; N];
        perm.copy_from_slice(&perm_vec);
        orbit.insert(pable.permute_by(&perm));
    }
    orbit
}
