use rustc_hash::{FxHashMap, FxHashSet};

use super::{
    dehydrogenation::generate_dehydrogenated,
    matrix::{
        AdjacencyBitMatrix, AdjacencyTwoBitsMatrix, Features, RedundantSaturatedHydrocarbonIter,
    },
    permutation::{calc_orbit_stabilizer, Permutation, PermutationStore},
};

pub fn create_feat2skeletons<const N: usize>() -> FxHashMap<Features<N>, Vec<AdjacencyBitMatrix<N>>>
{
    let mut feat2skeletons = FxHashMap::default();
    for (skeleton, feat) in RedundantSaturatedHydrocarbonIter::new() {
        feat2skeletons
            .entry(feat)
            .or_insert_with(Vec::new)
            .push(skeleton);
    }
    feat2skeletons
}

pub fn dehydrogenate_feat2skeletons<const N: usize>(
    feat2skeletons: FxHashMap<Features<N>, Vec<AdjacencyBitMatrix<N>>>,
    perm_store: &PermutationStore<N>,
) -> Vec<AdjacencyTwoBitsMatrix<N>> {
    let mut result = Vec::new();
    for (feat, skeletons) in feat2skeletons {
        let skeletons_stabilizers = remove_duplicates(skeletons.into_iter(), feat, perm_store);
        for (skeleton, stabilizer) in skeletons_stabilizers {
            let dehydrogenated = generate_dehydrogenated(skeleton.into(), &stabilizer);
            result.extend(dehydrogenated);
        }
    }
    result
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
