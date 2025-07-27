use std::collections::HashSet;

use super::{matrix::AdjacencyMatrix, permutation::calc_orbit};

pub fn generate_dehydrogenated<const N: usize>(
    skeleton: AdjacencyMatrix<N>,
) -> Vec<AdjacencyMatrix<N>> {
    let mut result = vec![skeleton];
    let mut queue = vec![skeleton];
    let mut seen_orbits = HashSet::new();

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
    hydrocarbon: &AdjacencyMatrix<N>,
) -> Vec<(usize, usize)> {
    let mut ables = [false; N];
    let mut pairs = Vec::with_capacity(2 * N);
    for irow in 0..N {
        if hydrocarbon.degree_of(irow) >= 4 {
            continue;
        }
        ables[irow] = true;

        // irow > icol の部分のみ検査
        for (icol, &able) in ables.iter().enumerate().take(irow) {
            if able && hydrocarbon.element_at(irow, icol) != 0 {
                pairs.push((irow, icol));
            }
        }
    }
    pairs
}
