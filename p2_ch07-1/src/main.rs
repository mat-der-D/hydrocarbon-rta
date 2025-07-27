#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod module;

use rustc_hash::FxHashMap;

use module::{
    create_feat2skeletons, dehydrogenate_feat2skeletons, AdjacencyTwoBitsMatrix, PermutationStore,
};

fn gen_all_hydrocarbons<const N: usize>() -> Vec<AdjacencyTwoBitsMatrix<N>> {
    let perm_store = PermutationStore::new();
    let feat2skeletons = create_feat2skeletons();
    dehydrogenate_feat2skeletons(feat2skeletons, &perm_store)
}

fn count_hydrogen<const N: usize>(hydrocarbon: &AdjacencyTwoBitsMatrix<N>) -> usize {
    let mut degrees_sum = 0;
    for i in 0..N {
        degrees_sum += hydrocarbon.degree_of(i) as usize;
    }
    4 * N - degrees_sum
}

fn run_impl<const N: usize>() {
    let hydrocarbons = gen_all_hydrocarbons::<N>();
    let mut counts = FxHashMap::default();
    for hydrocarbon in hydrocarbons {
        let num_h = count_hydrogen(&hydrocarbon);
        *counts.entry(num_h).or_insert(0) += 1;
    }
    println!("===== [C = {N:>2}] =====");
    println!("#H: #Hydrocarbons");
    for num_h in (0..=(N + 1)).map(|x| 2 * x) {
        println!("{:>2}: {}", num_h, counts.get(&num_h).unwrap_or(&0));
    }
}

macro_rules! run {
    ($($n:expr),+) => {
        $(run_impl::<$n>();)+
    };
}

fn main() {
    run!(2, 3, 4, 5, 6, 7, 8, 9, 10);
}
