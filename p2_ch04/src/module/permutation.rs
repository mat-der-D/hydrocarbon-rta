use std::hash::Hash;

use rustc_hash::FxHashSet;

pub trait Permutable<const N: usize> {
    fn permute_by(&self, perm: &[usize; N]) -> Self;
}

pub fn calc_orbit_stabilizer<const N: usize, P>(
    pable: P,
    perm_iter: impl Iterator<Item = [usize; N]>,
) -> (FxHashSet<P>, Vec<[usize; N]>)
where
    P: Permutable<N> + Eq + Hash,
{
    let mut orbit = FxHashSet::default();
    let mut stabilizer = Vec::new();
    for perm in perm_iter {
        let permuted = pable.permute_by(&perm);
        if permuted == pable {
            stabilizer.push(perm);
        }
        orbit.insert(permuted);
    }
    (orbit, stabilizer)
}

pub fn calc_orbit<const N: usize, P>(
    pable: P,
    perm_iter: impl Iterator<Item = [usize; N]>,
) -> FxHashSet<P>
where
    P: Permutable<N> + Eq + Hash,
{
    let mut orbit = FxHashSet::default();
    for perm in perm_iter {
        orbit.insert(pable.permute_by(&perm));
    }
    orbit
}
