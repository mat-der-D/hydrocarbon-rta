use std::hash::Hash;

use itertools::Itertools;
use rustc_hash::{FxHashMap, FxHashSet};

#[derive(Debug, Clone)]
pub struct PermutationStore<const N: usize> {
    memory: FxHashMap<u16, Vec<[usize; N]>>,
}

impl<const N: usize> PermutationStore<N> {
    pub fn new() -> Self {
        let mut memory = FxHashMap::default();
        for bits in (0..(1 << (N - 1))).map(|x| x << 1) {
            let key_array = Self::bits_to_key_array(bits);
            let perm_vec = Self::generate_perms(key_array);
            memory.insert(bits, perm_vec);
        }
        Self { memory }
    }

    pub fn get(&self, key: &u16) -> Option<&[[usize; N]]> {
        self.memory.get(key).map(|v| v.as_slice())
    }

    fn bits_to_key_array(bits: u16) -> [u8; N] {
        let mut key_array = [0; N];
        for (i, item) in key_array.iter_mut().enumerate() {
            *item = ((bits >> i) & 1) as u8;
        }
        key_array
    }

    fn generate_perms(key_array: [u8; N]) -> Vec<[usize; N]> {
        let mut perms_vec = Vec::new();
        Self::generate_perms_impl(&key_array, 0, &mut [0; N], &mut perms_vec);
        perms_vec
    }

    fn generate_perms_impl(
        key_array: &[u8],
        i_start: usize,
        perm: &mut [usize; N],
        perms_vec: &mut Vec<[usize; N]>,
    ) {
        if key_array.is_empty() {
            perms_vec.push(*perm);
            return;
        }

        let same_count = key_array.iter().take_while(|&&x| x == key_array[0]).count();
        let residual = &key_array[same_count..];
        let i_start_next = i_start + same_count;
        for partial_perm in (i_start..i_start_next).permutations(same_count) {
            perm[i_start..i_start_next].copy_from_slice(&partial_perm);
            Self::generate_perms_impl(residual, i_start_next, perm, perms_vec);
        }
    }
}

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
