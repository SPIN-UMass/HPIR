#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/LLL.h>
#include <NTL/vector.h>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <ctime>
#include <math.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <sys/time.h>
#include <NTL/mat_ZZ_p.h>
using namespace std;
NTL_CLIENT
void HPIR(){
  int size_of_prime_numbers= 512;
  uint32_t bytes_per_word = 512/8;
  uint32_t num_q;
  uint32_t words_per_block;
  uint32_t num_blocks;

  vec_ZZ primeNumbers;
  ZZ n=to_ZZ(1);
  primeNumbers.SetLength(2);
  for (size_t i = 0; i < 2; i++) {
    primeNumbers[i] = GenPrime_ZZ(size_of_prime_numbers,80);
    mul(n, n, primeNumbers[i]);
  }
  ZZ DBPrime;
  if (primeNumbers[0] < primeNumbers[1]){
    DBPrime = primeNumbers[0];
  }
  else{
    DBPrime = primeNumbers[1];
  }
  ZZ_p::init(n);
  ZZ mywrd_raw;



mat_ZZ_p tempDB;
double mult_rspp[] = {0.05, 0.1 , 0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0};
for (int iii=0; iii<9; iii++){

  int query_indeices[] = {1,2,3,4,5};
  ZZ prime0 =  primeNumbers[0];
  double total_time = 0;
  double server1_time = 0;
  double server2_time = 0;
  double client_time = 0;
  double Client_prepear_time = 0;
  double client_extract_time = 0;


  words_per_block = floor(sqrt((mult_rspp[iii]*1024*1024*1024)/bytes_per_word));//1536;//floor(sqrt(mult_rs));
  num_blocks = words_per_block;
  num_q = 3;
  tempDB.SetDims(num_blocks, words_per_block);

  for(uint32_t i=0; i<num_blocks; i++){
    for(uint32_t j=0; j<words_per_block; j++){
      tempDB[i][j] = to_ZZ_p(RandomBnd(DBPrime));
    }
  }

  ZZ_p myprod;
  ZZ_p mywrd;
  mat_ZZ_p tempmyResult;
  mat_ZZ_p myResult;
  mat_ZZ_p tempmyResult1;
  mat_ZZ_p myResult1;
  mat_ZZ_p tempmyResult2;
  mat_ZZ_p myResult2;
  mat_ZZ_p tempmyResult3;
  mat_ZZ_p myResult3;

  vec_ZZ_pX mypolys;
  mat_ZZ_p myQuery;
  mat_ZZ_p myQuery1;
  mat_ZZ_p myQuery2;
  clock_t t1 = clock();
  mat_ZZ_p A_y;
  ZZ_p temp_y;
  ZZ temprandomNumber_y;
  A_y.SetDims(num_blocks, num_q);
  clear(A_y);
  for (int i=0; i<num_blocks; i++){
    for (int j=0; j<num_q; j++){
      temprandomNumber_y = RandomBnd(n);
      conv(temp_y, temprandomNumber_y);
      if((temp_y != 0) && (temprandomNumber_y%primeNumbers[0]!=0) && (temprandomNumber_y%primeNumbers[1]!=0))
      A_y[i][j] = temp_y;
      else
      j--;
    }
  }

  // std::cout << "********************random point matrices******************" << '\n';
  mat_ZZ_p Points;
  ZZ_p temp_point;

  Points.SetDims(num_blocks, num_q);
  clear(Points);
  for (int i=0; i<num_blocks; i++){
    for (int j=0; j<num_q; j++){
      if(j==0){
        Points[i][j] = A_y[i][j];
      }
      else{
        if(j%2==0){
          conv(temp_point, primeNumbers[0]);
          mul(Points[i][j], A_y[i][j], temp_point);
        }else{
          conv(temp_point, primeNumbers[1]);
          mul(Points[i][j], A_y[i][j], temp_point);
        }
      }
    }
  }

  for (int i = 0; i < num_q-1; i++) {
    Points[query_indeices[i]][i+1] = Points[query_indeices[i]][i+1] + to_ZZ_p(1);//A_y[query_indeices[i]][i+1];
  }

  // std::cout << "*********************first Xcoordinates chosen randmly***************" << '\n';
  vec_ZZ_p Xcoordinates;
  ZZ_p temp2;
  ZZ temp2_ZZ;
  conv(temp2_ZZ, primeNumbers[0]);
  Xcoordinates.SetLength(num_q);
  clear(Xcoordinates);
  bool notRepeated = true;
  for(int i=0; i<num_q; i++){
    conv(temp2, RandomBnd(n));
    if(temp2 != 0){
      for (int j = 0; j < i; j++) {
        conv(temp2_ZZ, Xcoordinates[j]-temp2);
        if(temp2 == Xcoordinates[j])
        notRepeated = false;
        else{
          if(num_q>2){
            for(int m=0; m<2;m++){
              if(temp2_ZZ%primeNumbers[m]==0)
              notRepeated = false;
            }
          }
          else{
            for(int m=0; m<2;m++){
              if(temp2_ZZ%primeNumbers[m]==0)
              notRepeated = false;
            }
          }
        }
      }
      if(notRepeated){
        Xcoordinates[i] = temp2;
      }
      else{
        i--;
        notRepeated = true;
      }
    }
    else
    i--;
  }

  // std::cout << "*********************second Xcoordinates chosen randmly***************" << '\n';
  vec_ZZ_p Xcoordinates2;
  ZZ_p temp20;
  Xcoordinates2.SetLength(num_q);
  clear(Xcoordinates2);
  bool notRepeated2 = true;
  for(int i=0; i<num_q; i++){
    // temp2 = RandomBnd(prime0);
    conv(temp20, RandomBnd(n));
    if(temp20 != 0){
      for (int j = 0; j < i; j++) {
        conv(temp2_ZZ, Xcoordinates2[j]-temp20);
        if(temp20 == Xcoordinates2[j])
        notRepeated2 = false;
        else{
          if(num_q>2){
            for(int m=0; m<2;m++){
              if(temp2_ZZ%primeNumbers[m]==0)
              notRepeated2 = false;
            }
          }
          else{
            for(int m=0; m<2;m++){
              if(temp2_ZZ%primeNumbers[m]==0)
              notRepeated2 = false;
            }
          }
        }
      }
      for (int j = 0; j < num_q; j++) {
        if(temp20 == Xcoordinates[j])
        notRepeated2 = false;
      }
      if(notRepeated2){
        Xcoordinates2[i] = temp20;
      }
      else{
        i--;
        notRepeated2 = true;
      }
    }
    else
    i--;
  }

  // std::cout << "***********************real xcoordinates***************" << '\n';
  mat_ZZ_p final_x;
  ZZ_p temp_point_x;
  ZZ_p temp_point_xx;
  ZZ_p temp_point_xxx;
  ZZ_p temp_point_xxxx;
  ZZ temp23;
  bool checked_first_element = false;
  bool checked_second_element = false;

  final_x.SetDims(num_blocks, num_q);
  clear(final_x);
  for (int i=0; i<num_blocks; i++){
    for (int j=0; j<num_q; j++){
      if(checked_first_element){
        checked_first_element = false;
        j = 0;
      }
      if(j==0){
        conv(temp_point_xxx,RandomBnd(prime0));
        final_x[i][j] = temp_point_xxx;
        for(int k=0; k<num_q-1;k++){
          if(temp_point_xxx==Xcoordinates2[k]){
            checked_first_element = true;
          }
        }
      }
      else{
        if(j%2==0){
          conv(temp_point_xxx,RandomBnd(prime0));
          conv(temp_point_x, primeNumbers[0]);
          mul(temp_point_xx, temp_point_xxx, temp_point_x);
          add(temp_point_xxxx, temp_point_xx, Xcoordinates2[j-1]);
          for (size_t k = 0; k < j; k++) {
            conv(temp23, (temp_point_xxxx-final_x[i][k]));
            if(temp23%primeNumbers[0]==0 || temp23%primeNumbers[1]==0){
              checked_second_element = true;
            }
          }
          if(checked_second_element){
            checked_second_element = false;
            j--;
          }else{
            final_x[i][j] = temp_point_xxxx;
          }
        }else{
          conv(temp_point_xxx,RandomBnd(prime0));
          conv(temp_point_x, primeNumbers[1]);
          mul(temp_point_xx, temp_point_xxx, temp_point_x);
          add(temp_point_xxxx, temp_point_xx, Xcoordinates2[j-1]);
          for (int k = 0; k < j; k++) {
            conv(temp23, (temp_point_xxxx-final_x[i][k]));
            if(temp23%primeNumbers[0]==0 || temp23%primeNumbers[1]==0){
              checked_second_element = true;
            }
          }
          if(checked_second_element){
            checked_second_element = false;
            j--;
          }else{
            final_x[i][j] = temp_point_xxxx;
          }
        }
      }
    }
  }

  // std::cout << "*******************Generating functions******************************" << '\n';
  ZZ_p coff1;
  ZZ_p coff2;
  ZZ_p prime1_ZZP;
  mypolys.SetLength(num_blocks);
  vec_ZZ_p keys;
  keys.SetLength(num_q-1);

  for (int i=0; i<num_blocks; i++){
    interpolate(mypolys[i], final_x[i], Points[i]);
  }
  for (int j = 0; j < num_q-1; j++) {
    eval(keys[j], mypolys[query_indeices[j]], Xcoordinates2[j]);
  }

  myQuery1.SetDims(num_q-1, num_blocks);
  for (int j = 0; j < num_q-1; j++) {
    for (int i=0; i<num_blocks; i++){
      eval(myQuery1[j][i], mypolys[i], Xcoordinates[j]);
    }
  }

  myQuery2.SetDims(1, num_blocks);

  for (int i=0; i<num_blocks; i++){
    eval(myQuery2[0][i], mypolys[i], Xcoordinates[num_q-1]);
  }

  clock_t t5 = clock();

  // std::cout << "*************************serverside operations**************" << '\n';
  tempmyResult1.SetDims(num_q-1,words_per_block);
  myResult1.SetDims(words_per_block, num_q-1);

  tempmyResult1 = myQuery1 * tempDB;
  myResult1 = transpose(tempmyResult1);

  clock_t t55 = clock();

  tempmyResult2.SetDims(1,words_per_block);
  myResult2.SetDims(words_per_block, 1);

  tempmyResult2 = myQuery2 * tempDB;
  myResult2 = transpose(tempmyResult2);

  clock_t t6 = clock();
  tempmyResult.SetDims(num_q,words_per_block);
  myResult.SetDims(words_per_block, num_q);

  for (int i = 0; i < words_per_block; i++) {
    for (int j = 0; j <num_q-1 ; j++) {
      myResult[i][j] = myResult1[i][j];
    }
  }
  for (int i = 0; i < words_per_block; i++) {
    myResult[i][num_q-1] = myResult2[i][0];
  }

  // std::cout << "##################Interpolation#################" << '\n';
  vec_ZZ_pX phi;
  phi.SetLength(words_per_block);
  for (int i = 0; i < words_per_block; i++) {
    interpolate(phi[i], Xcoordinates, myResult[i]);
  }


  // std::cout << "##################Data Extraction#################" << '\n';
  ZZ_p temp5;
  ZZ temp8;
  ZZ temp6;
  ZZ temp7;
  unsigned char * mybytes = new unsigned char[bytes_per_word];
  for (size_t j = 0; j < num_q-1; j++) {
    for (size_t i = 0; i < words_per_block; i++) {
      if(j%2==1){
        eval(temp5, phi[i], Xcoordinates2[j]);
        conv(temp8,temp5);
        conv(temp7, keys[j]);
        InvMod(temp6, temp7%primeNumbers[0], primeNumbers[0]);
        temp7 = (temp8*temp6)%primeNumbers[0];
        BytesFromZZ(mybytes, temp7, bytes_per_word);
        conv(temp5, temp7);
        if (tempDB[query_indeices[j]][i] != temp5){
          std::cout << "$$$$$$$$$$$$$$$$$mismatch$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        }
      }else{
        eval(temp5, phi[i], Xcoordinates2[j]);
        conv(temp8,temp5);
        conv(temp7, keys[j]);
        InvMod(temp6, temp7%primeNumbers[1], primeNumbers[1]);
        temp7 = (temp8*temp6)%primeNumbers[1];
        BytesFromZZ(mybytes, temp7, bytes_per_word);
        conv(temp5, temp7);
        if (tempDB[query_indeices[j]][i] != temp5){
          std::cout << "$$$$$$$$$$$$$$$$$4mismatch$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        }
      }
    }
  }

  clock_t t7 = clock();
  total_time += double(t7-t6)+double(t55-t1);
  Client_prepear_time += double(t5-t1);
  client_extract_time += double(t7-t6);
  client_time += double((t7-t6)+(t5-t1));
  server1_time += double(t55-t5);
  server2_time += double(t6-t55);
  double database_size = double((num_blocks*words_per_block*bytes_per_word)/(1024));

  std::cout << "average total_time for: " << database_size <<" KB is: " << total_time << '\n';
  std::cout << "average server1_time for: " << database_size <<" KB is: " << server1_time << '\n';
  std::cout << "average server2_time for: " << database_size <<" KB is: " << server2_time << '\n';
  std::cout << "average Client_prepear_time for: " << database_size <<" KB is: " << Client_prepear_time << '\n';
  std::cout << "average client_extract_time for: " << database_size <<" KB is: " << client_extract_time << '\n';
  std::cout << "average client_time for: " << database_size <<" KB is: " << client_time << '\n';

  std::cout << "##############################finish##########" << '\n';
}
}
int main()
{
  HPIR();
  return 0;
}
