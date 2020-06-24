# HPIR
A repository for Heterogeneous Private Information Retrieval. HPIR allows a client to download an element from an untrusted servers in with heterogeneous communication/computation overhead. HPIR was introduced in [paper](https://people.cs.umass.edu/~amir/papers/NDSS20-HPIR.pdf) by [Hamid Mozaffari](https://github.com/lordhamiss) and Amir Houmansadr. 

# Compiling HPIR
HPIR depends on [NTL](https://www.shoup.net/ntl/). We have tested HPIR with NTL 11.4.3 and g++ 7.4.0.



`g++ -g -O2 -std=c++11 -pthread -march=native HPIR.cpp -o HPIR -lntl -lgmp -lm`
