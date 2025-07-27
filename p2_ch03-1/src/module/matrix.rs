use rustc_hash::FxHashSet;

use super::{
    hot_bit_iter::HotBitIter,
    permutation::{calc_orbit, Permutable},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AdjacencyBitMatrix<const N: usize> {
    rows: [u16; N],
}

impl<const N: usize> AdjacencyBitMatrix<N> {
    const ZERO: Self = Self { rows: [0; N] };

    fn element_at(&self, irow: usize, icol: usize) -> u16 {
        self.rows[irow] >> icol & 1
    }

    fn flip_at(&mut self, irow: usize, icol: usize) {
        self.rows[irow] ^= 1 << icol;
        self.rows[icol] ^= 1 << irow;
    }

    fn degree_of(&self, idx: usize) -> u32 {
        self.rows[idx].count_ones()
    }

    const REPUNIT: u16 = {
        let mut repunit = 0;
        let mut i = 0;
        while i < N {
            repunit |= 1 << i;
            i += 1;
        }
        repunit
    };

    fn is_connected(&self) -> bool {
        let mut visited = 1u16;
        for _ in 0..N {
            for irow in HotBitIter::from(visited) {
                visited |= self.rows[irow];
            }
            if visited == Self::REPUNIT {
                return true;
            }
        }
        false
    }
}

impl<const N: usize> Permutable<N> for AdjacencyBitMatrix<N> {
    fn permute_by(&self, perm: &[usize; N]) -> Self {
        let mut rows = [0; N];
        for (irow_old, &irow_new) in perm.iter().enumerate() {
            let row_old = self.rows[irow_old];
            let row_new = &mut rows[irow_new];

            for icol_old in HotBitIter::from(row_old) {
                let icol_new = perm[icol_old];
                *row_new |= 1 << icol_new;
            }
        }
        Self { rows }
    }
}

impl<const N: usize> std::fmt::Display for AdjacencyBitMatrix<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for irow in 0..N {
            for icol in 0..N {
                write!(f, "{} ", self.element_at(irow, icol))?;
            }
            if irow == N - 1 {
                break;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AdjacencyTwoBitsMatrix<const N: usize> {
    rows: [u32; N],
}

impl<const N: usize> AdjacencyTwoBitsMatrix<N> {
    pub fn element_at(&self, irow: usize, icol: usize) -> u32 {
        self.rows[irow] >> (2 * icol) & 0b11
    }

    pub fn increment_at(&mut self, irow: usize, icol: usize) {
        self.rows[irow] += 1 << (2 * icol);
        self.rows[icol] += 1 << (2 * irow);
    }

    pub fn degree_of(&self, idx: usize) -> u32 {
        let row = self.rows[idx];
        let odd_bits = row & 0x5555_5555; // 0b0101...0101
        let even_bits = row & 0xaaaa_aaaa; // 0b1010...1010
        odd_bits.count_ones() + 2 * even_bits.count_ones()
    }
}

impl<const N: usize> Permutable<N> for AdjacencyTwoBitsMatrix<N> {
    fn permute_by(&self, perm: &[usize; N]) -> Self {
        let mut rows = [0; N];
        for (irow_old, &irow_new) in perm.iter().enumerate() {
            let row_old = self.rows[irow_old];
            let row_new = &mut rows[irow_new];

            for num_zeros_old in HotBitIter::from(row_old) {
                let icol_old = num_zeros_old / 2;
                let icol_new = perm[icol_old];
                let num_zeros_new = icol_new * 2 + num_zeros_old % 2;
                *row_new |= 1 << num_zeros_new;
            }
        }
        Self { rows }
    }
}

impl<const N: usize> From<AdjacencyBitMatrix<N>> for AdjacencyTwoBitsMatrix<N> {
    fn from(mat: AdjacencyBitMatrix<N>) -> Self {
        let mut rows = [0; N];
        for (row_2, row_1) in rows.iter_mut().zip(mat.rows.into_iter()) {
            for icol in HotBitIter::from(row_1) {
                *row_2 |= 1 << (2 * icol);
            }
        }
        Self { rows }
    }
}

impl<const N: usize> std::fmt::Display for AdjacencyTwoBitsMatrix<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for irow in 0..N {
            for icol in 0..N {
                write!(f, "{} ", self.element_at(irow, icol))?;
            }
            if irow == N - 1 {
                break;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
struct Cursor<const N: usize> {
    irow: usize,
    icol: usize,
}

impl<const N: usize> Cursor<N> {
    fn new() -> Self {
        Self { irow: 0, icol: 1 }
    }

    const FIRST: (usize, usize) = (0, 1);
    const LAST: (usize, usize) = (N - 2, N - 1);

    fn is_at_first(&self) -> bool {
        (self.irow, self.icol) == Self::FIRST
    }

    fn is_at_last(&self) -> bool {
        (self.irow, self.icol) == Self::LAST
    }

    fn move_next(&mut self) {
        if self.is_at_last() {
            return;
        }

        (self.irow, self.icol) = if self.icol == N - 1 {
            (self.irow + 1, self.irow + 2)
        } else {
            (self.irow, self.icol + 1)
        }
    }

    fn move_prev(&mut self) {
        if self.is_at_first() {
            return;
        }

        (self.irow, self.icol) = if self.icol == self.irow + 1 {
            (self.irow - 1, N - 1)
        } else {
            (self.irow, self.icol - 1)
        }
    }
}

struct AdjacencyBitMatrixIter<const N: usize> {
    current: AdjacencyBitMatrix<N>,
    cursor: Cursor<N>,
}

impl<const N: usize> AdjacencyBitMatrixIter<N> {
    fn new() -> Self {
        Self {
            current: AdjacencyBitMatrix::ZERO,
            cursor: Cursor::new(),
        }
    }

    fn element_at_cursor(&self) -> u16 {
        self.current.element_at(self.cursor.irow, self.cursor.icol)
    }

    fn flip_at_cursor(&mut self) {
        self.current.flip_at(self.cursor.irow, self.cursor.icol);
    }
}

impl<const N: usize> Iterator for AdjacencyBitMatrixIter<N> {
    type Item = AdjacencyBitMatrix<N>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut forward = !self.cursor.is_at_last();
        loop {
            if forward {
                if self.cursor.is_at_last() {
                    return Some(self.current);
                }
                self.cursor.move_next();
            } else {
                self.flip_at_cursor();
                if self.element_at_cursor() == 0 {
                    if self.cursor.is_at_first() {
                        return None;
                    }
                    self.cursor.move_prev();
                } else {
                    forward = true;
                }
            }
        }
    }
}

pub struct SaturatedHydrocarbonIter<const N: usize> {
    mat_iter: AdjacencyBitMatrixIter<N>,
    seen_orbits: FxHashSet<AdjacencyBitMatrix<N>>,
}

impl<const N: usize> SaturatedHydrocarbonIter<N> {
    pub fn new() -> Self {
        Self {
            mat_iter: AdjacencyBitMatrixIter::new(),
            seen_orbits: FxHashSet::default(),
        }
    }

    fn is_hydrocarbon(mat: &AdjacencyBitMatrix<N>) -> bool {
        (0..N).all(|i| mat.degree_of(i) <= 4) && mat.is_connected()
    }
}

impl<const N: usize> Iterator for SaturatedHydrocarbonIter<N> {
    type Item = AdjacencyBitMatrix<N>;

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
