use rustc_hash::FxHashSet;

use super::permutation::calc_orbit;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AdjacencyMatrix<const N: usize> {
    elements: [[u8; N]; N],
}

impl<const N: usize> AdjacencyMatrix<N> {
    const ZERO: Self = Self {
        elements: [[0; N]; N],
    };

    pub fn element_at(&self, irow: usize, icol: usize) -> u8 {
        self.elements[irow][icol]
    }

    pub fn increment_at(&mut self, irow: usize, icol: usize) {
        self.elements[irow][icol] += 1;
        self.elements[icol][irow] += 1;
    }

    fn set_zero_at(&mut self, irow: usize, icol: usize) {
        self.elements[irow][icol] = 0;
        self.elements[icol][irow] = 0;
    }

    fn is_connected(&self) -> bool {
        let mut visited = [false; N];
        visited[0] = true;
        for _ in 0..N {
            let last_visited = visited;
            for icol in 0..N {
                if !last_visited[icol] {
                    continue;
                }
                for irow in 0..N {
                    if self.elements[irow][icol] != 0 {
                        visited[irow] = true;
                    }
                }
            }
            if visited.iter().all(|&v| v) {
                return true;
            }
        }
        false
    }

    pub fn permute_indices(&self, perm: &[usize; N]) -> Self {
        let mut elements = [[0; N]; N];
        for (irow_old, &irow_new) in perm.iter().enumerate() {
            for (icol_old, &icol_new) in perm.iter().enumerate() {
                elements[irow_new][icol_new] = self.elements[irow_old][icol_old];
            }
        }
        Self { elements }
    }

    pub fn degree_of(&self, idx: usize) -> u8 {
        self.elements[idx].iter().sum()
    }
}

impl<const N: usize> std::fmt::Display for AdjacencyMatrix<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for irow in 0..N {
            for icol in 0..N {
                write!(f, "{} ", self.elements[irow][icol])?;
            }
            if irow == N - 1 {
                break;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

struct AdjacencyMatrixIter<const N: usize> {
    current: AdjacencyMatrix<N>,
}

impl<const N: usize> AdjacencyMatrixIter<N> {
    fn new() -> Self {
        Self {
            current: AdjacencyMatrix::ZERO,
        }
    }
}

impl<const N: usize> Iterator for AdjacencyMatrixIter<N> {
    type Item = AdjacencyMatrix<N>;

    fn next(&mut self) -> Option<Self::Item> {
        let (mut irow, mut icol) = (0, 1);
        loop {
            if self.current.element_at(irow, icol) == 0 {
                self.current.increment_at(irow, icol);
                return Some(self.current);
            } else {
                self.current.set_zero_at(irow, icol);
                if irow == N - 2 {
                    return None;
                }
            }

            (irow, icol) = if icol < N - 1 {
                (irow, icol + 1)
            } else {
                (irow + 1, irow + 2)
            };
        }
    }
}

pub struct SaturatedHydrocarbonIter<const N: usize> {
    mat_iter: AdjacencyMatrixIter<N>,
    seen_orbits: FxHashSet<AdjacencyMatrix<N>>,
}

impl<const N: usize> SaturatedHydrocarbonIter<N> {
    pub fn new() -> Self {
        Self {
            mat_iter: AdjacencyMatrixIter::new(),
            seen_orbits: FxHashSet::default(),
        }
    }

    fn is_hydrocarbon(mat: &AdjacencyMatrix<N>) -> bool {
        (0..N).all(|i| mat.degree_of(i) <= 4) && mat.is_connected()
    }
}

impl<const N: usize> Iterator for SaturatedHydrocarbonIter<N> {
    type Item = AdjacencyMatrix<N>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mat = self.mat_iter.next()?;
            if !Self::is_hydrocarbon(&mat) {
                continue;
            }
            if self.seen_orbits.contains(&mat) {
                continue;
            }
            self.seen_orbits.extend(calc_orbit(mat));
            return Some(mat);
        }
    }
}
