use rustc_hash::FxHashSet;

use super::{hot_bit_iter::HotBitIter, matrix::AdjacencyTwoBitsMatrix, permutation::calc_orbit};

pub fn generate_dehydrogenated<const N: usize>(
    skeleton: AdjacencyTwoBitsMatrix<N>,
) -> Vec<AdjacencyTwoBitsMatrix<N>> {
    let mut result = vec![skeleton];
    let mut queue = vec![skeleton];
    let mut seen_orbits = FxHashSet::default();

    while !queue.is_empty() {
        let mut next_queue = Vec::new();
        for hydrocarbon in queue {
            let index_pairs = find_possible_index_pairs(&hydrocarbon);
            for (irow, icol) in index_pairs {
                let mut dehydrogenated = hydrocarbon;
                dehydrogenated.increment_at(irow, icol);
                if seen_orbits.contains(&dehydrogenated) {
                    continue;
                }
                result.push(dehydrogenated);
                next_queue.push(dehydrogenated);
                seen_orbits.extend(calc_orbit(dehydrogenated));
            }
        }
        queue = next_queue;
        seen_orbits.clear();
    }
    result
}

fn find_possible_index_pairs<const N: usize>(
    hydrocarbon: &AdjacencyTwoBitsMatrix<N>,
) -> Vec<(usize, usize)> {
    let mut ables = 0u16;
    let mut pairs = Vec::with_capacity(2 * N);
    let max_degree = if N == 2 { 3 } else { 4 }; // N = 2 でのオーバーフロー防止
    for irow in 0..N {
        if hydrocarbon.degree_of(irow) >= max_degree {
            continue;
        }
        ables |= 1 << irow;

        for icol in HotBitIter::from(ables) {
            if hydrocarbon.element_at(irow, icol) != 0 {
                pairs.push((irow, icol));
            }
        }
    }
    pairs
}
