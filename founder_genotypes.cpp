#include <Rcpp.h>
#include <random>
#include <fstream>
#include <string>

//going to save the founder genos directly to file
//do the shuffling in the convert function
//now deciding how to do the poisson sampling

// the poolsize is the initial pool of individual haplotypes to be simulated
// afterwards the founder genotypes are obtained by splicing together the genotypes from this pool
//      so a smaller pool will results in a more related sample

// lambda is the lambda param for poisson distribution. 
//      will select a haplotype from the pool of n consecutive markers sampled from poisson distribution
//      larger lambda will result in longer haplotypes (more relatedness)

/// [[Rcpp::export]]
// void founder_genos_related(int poolsize, Rcpp::IntegerVector founders, int n_marker, double MAF, int lambda, std::string out) {
//     const int n_founders = founders.size(); 
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::poisson_distribution<> pois(lambda);
//     std::uniform_int_distribution<> unif(0, poolsize);

// 	std::ofstream outfile(out + ".ped");

//     //this vector contains all the genotypes for the n individuals in the initial pool
//     //stored as poolsize runs of size n_marker
//     //std vector default initializes to false
//     std::vector<char> pool_geno(poolsize * n_marker, '0'); 
//     std::fill_n(pool_geno.begin(), (int)(pool_geno.size() * MAF), '1');
//     std::shuffle(pool_geno.begin(), pool_geno.end(), gen);

//     bool done = false;
//     int position =0, run_len, sample_pool;

//     for (int i=0; i<n_founders; i++){
//         //position is number of markers
//         Rcout << i << "\n" ;
//         std::vector<char> founder_haplo1( n_marker);
//         std::vector<char> founder_haplo2( n_marker);
//         done = false;
//         position = 0;

//         while (not done){
//             run_len  = pois(gen);
//             sample_pool = unif(gen) * poolsize + position;
//             if (position +run_len > n_marker){
//                 done = true;
//                 run_len = n_marker - position;}
//             std::copy(pool_geno.begin() + sample_pool, pool_geno.begin() + sample_pool + run_len, founder_haplo1.begin() + position );
//             position = position + run_len;
//         }
//         done = false;
//         position = 0;
//         while (not done){
//             run_len  = pois(gen);
//             sample_pool = unif(gen) * poolsize + position;
//             if (position +run_len > n_marker){
//                 done = true;
//                 run_len = n_marker - position;}
//             std::copy(pool_geno.begin() + sample_pool, pool_geno.begin() + sample_pool + run_len, founder_haplo2.begin() + position );
//             position = position + run_len;
//         }  
        
//         //put the first 6 col of pedfile
//         outfile << "0\t" << founders[i] << "\t0\t0\t-9\t-9\t";
//         for (int j = 0; j<n_marker; j++){
//             outfile << founder_haplo1[j] << '\t' << founder_haplo2[j] << '\t';
//         }
//         outfile << std::endl;
//     }	
// }

// [[Rcpp::export]]
void founder_genos(Rcpp::IntegerVector founders, Rcpp::NumericVector MAF, std::string outfile) {
    const int n_founders = founders.size(); 
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::poisson_distribution<> pois(lambda);
    std::uniform_int_distribution<> unif(0, poolsize);

	std::ofstream outfile(out + ".ped");

    //this vector contains all the genotypes for the n individuals in the initial pool
    //stored as poolsize runs of size n_marker
    //std vector default initializes to false
    std::vector<char> pool_geno(poolsize * n_marker, '0'); 
    std::fill_n(pool_geno.begin(), (int)(pool_geno.size() * MAF), '1');
    std::shuffle(pool_geno.begin(), pool_geno.end(), gen);

    bool done = false;
    int position =0, run_len, sample_pool;

    for (int i=0; i<n_founders; i++){
        //position is number of markers
        Rcout << i << "\n" ;
        std::vector<char> founder_haplo1( n_marker);
        std::vector<char> founder_haplo2( n_marker);
        done = false;
        position = 0;

        while (not done){
            run_len  = pois(gen);
            sample_pool = unif(gen) * poolsize + position;
            if (position +run_len > n_marker){
                done = true;
                run_len = n_marker - position;}
            std::copy(pool_geno.begin() + sample_pool, pool_geno.begin() + sample_pool + run_len, founder_haplo1.begin() + position );
            position = position + run_len;
        }
        done = false;
        position = 0;
        while (not done){
            run_len  = pois(gen);
            sample_pool = unif(gen) * poolsize + position;
            if (position +run_len > n_marker){
                done = true;
                run_len = n_marker - position;}
            std::copy(pool_geno.begin() + sample_pool, pool_geno.begin() + sample_pool + run_len, founder_haplo2.begin() + position );
            position = position + run_len;
        }  
        
        //put the first 6 col of pedfile
        outfile << "0\t" << founders[i] << "\t0\t0\t-9\t-9\t";
        for (int j = 0; j<n_marker; j++){
            outfile << founder_haplo1[j] << '\t' << founder_haplo2[j] << '\t';
        }
        outfile << std::endl;
    }	
}