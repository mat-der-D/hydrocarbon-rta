use std::{collections::VecDeque, hash::Hash};

use rustc_hash::{FxHashMap, FxHashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Permutation<const N: usize> {
    raw: [usize; N],
}

impl<const N: usize> Permutation<N> {
    const fn new(raw: [usize; N]) -> Self {
        Self { raw }
    }

    const INDEX_ARRAY: [usize; N] = {
        let mut array = [0; N];
        let mut i = 0;
        while i < N {
            array[i] = i;
            i += 1;
        }
        array
    };

    const IDENTITY: Self = Self::new(Self::INDEX_ARRAY);

    fn new_cyclic(start: usize, end: usize) -> Self {
        let mut raw = Self::INDEX_ARRAY;
        for i in start..end {
            raw.swap(i, i + 1);
        }
        Self::new(raw)
    }

    fn inverse(&self) -> Self {
        let mut inverse = [0; N];
        for i in 0..N {
            inverse[self.raw[i]] = i;
        }
        Self::new(inverse)
    }

    fn permute<P: Permutable<N>>(&self, permutable: &P) -> P {
        permutable.permute_by(&self.raw)
    }
}

impl<const N: usize> std::ops::Mul for Permutation<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut raw = [0; N];
        for i in 0..N {
            raw[i] = self.raw[rhs.raw[i]];
        }
        Self::new(raw)
    }
}

#[derive(Debug, Clone)]
pub struct PermutationStore<const N: usize> {
    memory: FxHashMap<u16, Vec<Permutation<N>>>,
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

    pub fn get(&self, key: &u16) -> Option<&[Permutation<N>]> {
        self.memory.get(key).map(|v| v.as_slice())
    }

    fn bits_to_key_array(bits: u16) -> [u8; N] {
        let mut key_array = [0; N];
        for (i, item) in key_array.iter_mut().enumerate() {
            *item = ((bits >> i) & 1) as u8;
        }
        key_array
    }

    fn generate_perms(key_array: [u8; N]) -> Vec<Permutation<N>> {
        let mut perms = Vec::new();
        let mut i_start = 0;
        for i_end in 0..N {
            if i_end != N - 1 && key_array[i_start] == key_array[i_end + 1] {
                continue;
            }
            match i_end - i_start {
                0 => (),
                1 => perms.push(Permutation::new_cyclic(i_start, i_end)),
                2.. => {
                    perms.push(Permutation::new_cyclic(i_start, i_start + 1));
                    perms.push(Permutation::new_cyclic(i_start, i_end));
                }
            }
            i_start = i_end + 1;
        }
        perms
    }
}

pub trait Permutable<const N: usize> {
    fn permute_by(&self, perm: &[usize; N]) -> Self;
}

pub fn calc_orbit_stabilizer<const N: usize, P>(
    pable: P,
    generators: &[Permutation<N>],
    orbit_pre_alloc: usize,
) -> (impl IntoIterator<Item = P>, Vec<Permutation<N>>)
where
    P: Permutable<N> + Clone + Eq + Hash,
{
    let mut stabilizer = FxHashSet::default();

    let mut history = FxHashMap::with_capacity_and_hasher(orbit_pre_alloc, Default::default());
    history.insert(pable.clone(), Permutation::IDENTITY);

    let mut queue = VecDeque::new();
    queue.push_back((pable, Permutation::IDENTITY));

    while let Some((pable_tmp, p)) = queue.pop_front() {
        for &g in generators {
            let permuted = g.permute(&pable_tmp);
            let gp = g * p;
            if let Some(&p_hist) = history.get(&permuted) {
                if gp != p_hist {
                    stabilizer.insert(p_hist.inverse() * gp);
                }
            } else {
                history.insert(permuted.clone(), gp);
                queue.push_back((permuted, gp));
            }
        }
    }

    (history.into_keys(), stabilizer.into_iter().collect())
}

pub fn calc_orbit<const N: usize, P>(
    pable: P,
    generators: &[Permutation<N>],
) -> impl IntoIterator<Item = P>
where
    P: Permutable<N> + Clone + Eq + Hash,
{
    let mut orbit = FxHashSet::default();
    orbit.insert(pable.clone());

    let mut queue = VecDeque::new();
    queue.push_back(pable);

    while let Some(pable_tmp) = queue.pop_front() {
        for g in generators {
            let permuted = g.permute(&pable_tmp);
            if !orbit.contains(&permuted) {
                orbit.insert(permuted.clone());
                queue.push_back(permuted);
            }
        }
    }
    orbit
}
