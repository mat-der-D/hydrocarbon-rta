use std::{sync::mpsc, thread};

use rustc_hash::{FxHashMap, FxHashSet};

use super::{
    dehydrogenation::generate_dehydrogenated,
    matrix::{
        AdjacencyBitMatrix, AdjacencyTwoBitsMatrix, Features, RedundantSaturatedHydrocarbonIter,
    },
    permutation::{calc_orbit_stabilizer, Permutation, PermutationStore},
};

pub fn create_feat2skeletons<const N: usize>(
    fixed_digits: usize,
) -> FxHashMap<Features<N>, Vec<AdjacencyBitMatrix<N>>> {
    let skeleton_iters = RedundantSaturatedHydrocarbonIter::new_with_fixed_digits(fixed_digits);
    let (sender, receiver) = mpsc::channel();
    let mut result = FxHashMap::default();

    thread::scope(|s| {
        for iter in skeleton_iters {
            let sender = mpsc::Sender::clone(&sender);
            s.spawn(move || {
                for x in iter {
                    sender.send(x).unwrap();
                }
            });
        }
        drop(sender); // 最初の1個が余るので手動で drop

        for (skeleton, feat) in receiver {
            result.entry(feat).or_insert_with(Vec::new).push(skeleton);
        }
    });

    result
}

pub fn dehydrogenate_feat2skeletons<const N: usize>(
    feat2skeletons: FxHashMap<Features<N>, Vec<AdjacencyBitMatrix<N>>>,
    perm_store: &PermutationStore<N>,
    max_num_feats: usize,
) -> Vec<AdjacencyTwoBitsMatrix<N>> {
    let f2s = &feat2skeletons;
    let num_threads = f2s.len().div_ceil(max_num_feats);

    let mut hydrocarbons = Vec::new();
    let (sender, receiver) = mpsc::channel();

    thread::scope(|s| {
        for ith in 0..num_threads {
            let sender = mpsc::Sender::clone(&sender);
            s.spawn(move || {
                for (&feat, skeletons) in f2s.iter().skip(ith).step_by(num_threads) {
                    let skeletons_stabilizers =
                        remove_duplicates(skeletons.iter().copied(), feat, perm_store);
                    for (skeleton, stabilizer) in skeletons_stabilizers {
                        let dehydrogenated = generate_dehydrogenated(skeleton.into(), &stabilizer);
                        sender.send(dehydrogenated).unwrap();
                    }
                }
            });
        }
        drop(sender); // 最初の1個が余るので手動で drop

        for dehydrogenated in receiver {
            hydrocarbons.extend(dehydrogenated);
        }
    });

    hydrocarbons
}

fn remove_duplicates<const N: usize>(
    skeletons: impl Iterator<Item = AdjacencyBitMatrix<N>>,
    feat: Features<N>,
    perm_store: &PermutationStore<N>,
) -> Vec<(AdjacencyBitMatrix<N>, Vec<Permutation<N>>)> {
    let mut skeletons_stabilizers = Vec::new();
    let mut seen_orbits = FxHashSet::default();
    let generators = perm_store.get(&feat.make_key()).unwrap();
    for skeleton in skeletons {
        if seen_orbits.contains(&skeleton) {
            continue;
        }
        let (orbit, stabilizer) = calc_orbit_stabilizer(skeleton, generators);
        seen_orbits.extend(orbit);
        skeletons_stabilizers.push((skeleton, stabilizer));
    }
    skeletons_stabilizers
}
