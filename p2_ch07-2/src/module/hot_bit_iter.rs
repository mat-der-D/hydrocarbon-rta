#[derive(Debug)]
pub struct HotBitIter<T> {
    bits: T,
}

macro_rules! impl_hot_bit_iter {
    ($($t:ty),*) => {
        $(
            impl From<$t> for HotBitIter<$t> {
                fn from(bits: $t) -> Self {
                    Self { bits }
                }
            }

            impl Iterator for HotBitIter<$t> {
                type Item = usize;

                fn next(&mut self) -> Option<Self::Item> {
                    if self.bits == 0 {
                        return None;
                    }
                    let trailing_zeros = self.bits.trailing_zeros() as usize;
                    self.bits &= self.bits - 1;
                    Some(trailing_zeros)
                }
            }
        )*
    };
}

impl_hot_bit_iter!(u16, u32);
