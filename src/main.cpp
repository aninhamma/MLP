#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <sys/timeb.h>
#include <sys/resource.h>
#include <ctime>

using namespace std;

double **matrizAdj; // matriz de adjacencia
int dimension;      // quantidade total de vertices

void printData();

vector<int> solucaoFinal;

double tempo_construction = 0;
double tempo_swap = 0;
double tempo_reinsertion = 0;
double tempo_2opt = 0;
double tempo_orOpt2 = 0;
double tempo_orOpt3 = 0;
double tempo_perturbacao = 0;

double cpuTime()
{
  static struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return ((double)usage.ru_utime.tv_sec) + (((double)usage.ru_utime.tv_usec) / ((double)1000000));
}

void printTime()
{
  cout << "Tempo medio de execucao da SI: " << (tempo_construction / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao da troca: " << (tempo_swap / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao do Or-opt: " << (tempo_reinsertion / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao do Or-opt2: " << (tempo_orOpt2 / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao do Or-opt3: " << (tempo_orOpt3 / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao do 2-opt: " << (tempo_2opt / 10) << " (s)";
  cout << "\n"
       << "Tempo medio de execucao da perturbacao: " << (tempo_perturbacao / 10) << " (s)";

  cout << "\n\n";
}

struct InsertionInfo
{
  int noInserido;     // k a ser inserido
  int arestaRemovida; // aresta {i,j} na qual o k vai ser inserido
  double custo;       // delta
};

struct Subsequence
{
  double T, C;     // duracao e custo acumulado
  int W;           // custo de atraso
  int first, last; // primeiro e ulimo no da subsequencia

  inline static Subsequence Concatenate(Subsequence &sigma_1, Subsequence &sigma_2)
  {
    Subsequence sigma;
    double temp = matrizAdj[sigma_1.last][sigma_2.first];
    sigma.W = sigma_1.W + sigma_2.W;
    sigma.T = sigma_1.T + temp + sigma_2.T;
    sigma.C = sigma_1.C + sigma_2.W * (sigma_1.T + temp) + sigma_2.C;
    sigma.first = sigma_1.first;
    sigma.last = sigma_2.last;

    return sigma;
  }
};

void printSolucao(vector<int> sequencia)
{
  for (int i = 0; i < sequencia.size() - 1; i++)
    cout << sequencia[i] << "->";
  cout << sequencia.back() << endl;
}

int random(int low, int high)
{
  int valorRetorno = low + rand() % (high - low + 1);
  // cout << "low: " << low << endl;
  // cout << "high: " << high << endl;

  return valorRetorno;
}

void removeDaLista(vector<int> &CL, int n)
{ // funcao que remove o no inserido da lista de candidatos

  CL.erase(CL.begin() + n);
}

void escolher3NosAleatorios(vector<int> &sequencia, vector<int> &CL)
{ // funcao que escolhe os 3 primeiros nos da solucao(uma sequencia inicial)

  int n = 0; // comecando o no com o primeiro elemento inserido na sequencia
  int i;

  sequencia = {1, 1}; // inicilaizei a sequencia com a primeira cidade

  removeDaLista(CL, n); // remove o primeiro no inserido

  for (i = 0; i < 3; i++)
  { // for para repetir o processo de escolha aleatoria para pegar tres nos

    n = rand() % CL.size(); // escolhe aleatoriamente um no da lista de candidatos

    sequencia.insert(sequencia.end() - 1, CL[n]); // insere o no escolhido na penultima posicao

    removeDaLista(CL, n); // remove o no inserido da lista de candidados
  }
}

vector<InsertionInfo> calcularCustoInsercao(vector<int> &sequencia, vector<int> &CL)
{

  vector<InsertionInfo> custoInsercao((sequencia.size() - 1) * CL.size());

  int l = 0;
  for (int a = 0, b = 1; a < sequencia.size() - 1; a++, b++)
  {
    int i = sequencia[a];
    int j = sequencia[b];

    for (auto k : CL)
    {
      // cout << i << j << k << l;
      custoInsercao[l].custo = matrizAdj[i][k] + matrizAdj[i][j] - matrizAdj[i][j];
      custoInsercao[l].noInserido = k;
      custoInsercao[l].arestaRemovida = a;
      l++;
    }
  }
  return custoInsercao;
}

bool compares(InsertionInfo a, InsertionInfo b)
{
  return a.custo < b.custo;
}

void ordenarEmOrdemCrescente(vector<InsertionInfo> &custoInsercao)
{
  sort(custoInsercao.begin(), custoInsercao.end(), compares);
}

void InserirNaSolucao(vector<int> &sequencia, int selecionado, int arestaRemovida)
{

  sequencia.insert(sequencia.begin() + arestaRemovida + 1, selecionado);
}

double custoDaSolucao(vector<int> &sequencia)
{
  double custo = 0;
  for (int i = 0, j = 1; i < sequencia.size() - 1; i++, j++)
  {
    custo += matrizAdj[sequencia[i]][sequencia[j]];
  }
  return custo;
}

vector<int> construcao(double alpha)
{

  double inicioconstruction = cpuTime();
  vector<int> CL;
  vector<int> sequencia; // primeira solucao
  // vector <int> best;//solucao melhorada
  int i;
  double custo = 0;

  for (i = 0; i < dimension; i++)
  {
    CL.push_back(i + 1);
  }

  // cout << "Lista de candidatos: " << endl;

  // printSolucao(CL);

  escolher3NosAleatorios(sequencia, CL);

  // cout << "Tres cidades iniciais: " << endl;

  // printSolucao(sequencia);

  // cout << alpha << endl;

  vector<InsertionInfo> custoInsercao = calcularCustoInsercao(sequencia, CL);

  while (!CL.empty())
  {

    ordenarEmOrdemCrescente(custoInsercao);

    int custoS = (int)(ceil(alpha * custoInsercao.size()));

    // cout << custoInsercao.size() << endl;

    // cout << custoS << endl;

    int selecionado = rand() % custoS;

    // cout << selecionado << endl;

    InserirNaSolucao(sequencia, custoInsercao[selecionado].noInserido, custoInsercao[selecionado].arestaRemovida);

    for (int i = 0; i < CL.size(); i++)
    {
      int k = custoInsercao[selecionado].noInserido;
      if (CL[i] == k)
      {
        removeDaLista(CL, i);
      }
    }

    custoInsercao = calcularCustoInsercao(sequencia, CL);
  }

  // cout << "Solucao inicial: " << endl;

  // printSolucao(sequencia);

  double fimConstruction = cpuTime();
  tempo_construction += (fimConstruction - inicioconstruction);

  custo = custoDaSolucao(sequencia);

  // cout << "custo da solucao inicial: " << custo << endl;

  return sequencia;
}

void UpdateAllSubseq(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{

  int n = sequencia.size();

  // subsequencias de um unico no
  for (int i = 0; i < n; i++)
  {
    int v = sequencia[i];
    subseq_matrix[i][i].W = (i > 0);
    subseq_matrix[i][i].C = 0;
    subseq_matrix[i][i].T = 0;
    subseq_matrix[i][i].first = sequencia[i];
    subseq_matrix[i][i].last = sequencia[i];
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j - 1], subseq_matrix[j][j]);
    }
  }

  for (int i = n - 1; i >= 0; i--)
  {
    for (int j = i - 1; j >= 0; j--)
    {
      subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j + 1], subseq_matrix[j][j]);
    }
  }
}

/*inline double calculateSwapCost(int i, int j, vector <int> &sequencia){ //i e j sao cidades que serao trocadas

  double delta;

  double parteDaSeqInicial = -(matrizAdj[sequencia[i]][sequencia[i + 1]] + matrizAdj[sequencia[i]][sequencia[i - 1]] +
                              matrizAdj[sequencia[j]][sequencia[j + 1]] + matrizAdj[sequencia[j]][sequencia[j - 1]]);//parte da sequencia inicial que sera retirada quando as cidades trocadas nao forem vizinhas

  if(j == i + 1){ //para quando as cidades trocadas forem vizinhas
    delta = -matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[j]][sequencia[j + 1]] +
            matrizAdj[sequencia[i - 1]][sequencia[j]] + matrizAdj[sequencia[i]][sequencia[j + 1]];
  }else{
    delta = parteDaSeqInicial + matrizAdj[sequencia[i]][sequencia[j + 1]] + matrizAdj[sequencia[i]][sequencia[j - 1]] + matrizAdj[sequencia[j]][sequencia[i + 1]] + matrizAdj[sequencia[j]][sequencia[i - 1]];
  }

  return delta;

}*/

bool bestImprovementSwap(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicioSwap = cpuTime();
  double melhorCusto = subseq_matrix[0][dimension].C;
  int tamanho = dimension + 1;
  Subsequence sigma_1, sigma_2, sigma_3, last_sigma;
  // vector<int> copia = sequencia;

  int best_i = 0, best_j = 0;
  for (int i = 1; i < tamanho - 1; i++)
  {
    for (int j = i + 1; j < tamanho - 1; j++)
    {

      if (j == i + 1)
      { // para quando as cidades trocadas forem vizinhas
        sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i - 1], subseq_matrix[j][j]);
        sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[i][i]);
        last_sigma = Subsequence::Concatenate(sigma_2, subseq_matrix[j + 1][dimension]);
      }
      else{
        sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i - 1], subseq_matrix[j][j]);
        sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[i+1][j-1]);
        sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[i][i]);
        last_sigma = Subsequence::Concatenate(sigma_3, subseq_matrix[i+1][dimension]);
      }

      // delta = calculateSwapCost(i, j, sequencia);
      // custoTeste = custoDaSolucao(copia);

      /*if(custo + delta != custoTeste){
        cout << custo + delta << endl;
        cout << custoTeste << endl;
        cout << "i: " << i << "j: " << j << endl;
      }*/

      // copia = sequencia;

      if (last_sigma.C < melhorCusto)
      {
        // cout << "comparacao dos deltas" << endl;
        melhorCusto = last_sigma.C;
        best_i = i;
        best_j = j;
      }
    }
  }

  if (melhorCusto < subseq_matrix[0][dimension].C)
  {
    swap(sequencia[best_i], sequencia[best_j]);
    UpdateAllSubseq(sequencia, subseq_matrix);
    // cout << "custo antes da troca: " << custo << endl;
    //custo = custo + bestDelta;
    // cout << "custo apos a troca: " << custo << endl;
    double fimSwap = cpuTime();
    tempo_swap += (fimSwap - inicioSwap);
    return true;
  }
  else
  {
    double fimSwap = cpuTime();
    tempo_swap += (fimSwap - inicioSwap);
    return false;
  }
}

/*inline double calculate2OptCost(int i, int j, vector<int> &sequencia){
  double delta;

  delta = -matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[j]][sequencia[j + 1]] + matrizAdj[sequencia[i - 1]][sequencia[j]] + matrizAdj[sequencia[i]][sequencia[j + 1]];

  return delta;
}*/

bool bestImprovement2Opt(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicio2opt = cpuTime();
  double melhorCusto = subseq_matrix[0][dimension].C;
  int tamanho = dimension + 1;
  Subsequence sigma_1, sigma_2;
  // vector<int> copia = sequencia;

  int best_i = 0, best_j = 0;
  for (int i = 1; i < tamanho - 1; i++){
    for (int j = i + 1; j < tamanho - 1; j++){
      // reverse(copia.begin() + i, copia.begin() + j + 1);
      sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][i]);
      sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[j+1][dimension]);

      // delta = calculate2OptCost(i, j, sequencia);
      // custoTeste = custoDaSolucao(copia);

      /*if(custo + delta != custoTeste){
        cout << custo + delta << endl;
        cout << custoTeste << endl;
        cout << "i: " << i << "j: " << j << endl;
      }*/

      // opia = sequencia;

      if (sigma_2.C < melhorCusto)
      {
        melhorCusto = sigma_2.C;
        best_i = i;
        best_j = j;
      }
    }
  }

  if (melhorCusto < subseq_matrix[0][dimension].C)
  {
    reverse(sequencia.begin() + best_i, sequencia.begin() + best_j + 1);
    UpdateAllSubseq(sequencia, subseq_matrix);
    // cout << "custo antes do 2opt: " << custo << endl;
    //custo = custo + bestDelta;
    // cout << "custo apos o 2opt: " << custo << endl;
    double fim2opt = cpuTime();
    tempo_2opt += (fim2opt - inicio2opt);
    return true;
  }
  else
  {
    double fim2opt = cpuTime();
    tempo_2opt += (fim2opt - inicio2opt);
    return false;
  }
}

/*inline double calculateReinsertionCost(int i, int j, vector<int> &sequencia){
  double delta;

  if(i < j){
    delta = -matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i]][sequencia[i + 1]] - matrizAdj[sequencia[j]][sequencia[j + 1]] + matrizAdj[sequencia[i]][sequencia[j]] + matrizAdj[sequencia[i]][sequencia[j + 1]] + matrizAdj[sequencia[i - 1]][sequencia[i + 1]];
  }else{
    delta = -matrizAdj[sequencia[j]][sequencia[j - 1]] - matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i]][sequencia[i + 1]] + matrizAdj[sequencia[i]][sequencia[j]] + matrizAdj[sequencia[i - 1]][sequencia[i + 1]] + matrizAdj[sequencia[j - 1]][sequencia[i]];
  }
  return delta;
}*/

bool bestImprovementReinsertion(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicioreinsertion = cpuTime();
  double melhorCusto = subseq_matrix[0][dimension].C;
  //vector<int> copia = sequencia;
  int tamanho = dimension + 1;
  Subsequence sigma_1, sigma_2, sigma_3;

  int best_i = 0, best_j = 0;
  for (int i = 1; i < tamanho - 1; i++){
    for (int j = 1; j < tamanho - 1; j++){

      if (i != j){
        /*copia.erase(copia.begin() + i);
        copia.insert(copia.begin() + j, sequencia[i]);*/
        if (i < j)
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+1][j]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[i][i]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[j+1][dimension-1]);  
        }
        else
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][j-1], subseq_matrix[i][i]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[j][i-1]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[i+1][dimension-1]);  
        }
        // delta = calculateReinsertionCost(i, j, sequencia);
        // custoTeste = custoDaSolucao(copia);

        /*if(custo + delta != custoTeste){
          cout << custo + delta << endl;
          cout << custoTeste << endl;
          cout << "i: " << i << "j: " << j << endl;
        }*/

        // copia = sequencia;

        if (sigma_3.C < melhorCusto)
        {
          melhorCusto = sigma_3.C;
          best_i = i;
          best_j = j;
        }
      }
    }
  }

  if (melhorCusto < subseq_matrix[0][dimension].C)
  {
    int reinserir = sequencia[best_i];
    sequencia.erase(sequencia.begin() + best_i);
    sequencia.insert(sequencia.begin() + best_j, reinserir);
    UpdateAllSubseq(sequencia, subseq_matrix);

    // cout << "custo antes do reinsertion: " << custo << endl;
    //custo = custo + bestDelta;
    // cout << "custo apos o reinsertion: " << custo << endl;
    double fimReinsertion = cpuTime();
    tempo_reinsertion += (fimReinsertion - inicioreinsertion);
    return true;
  }
  else
  {
    double fimReinsertion = cpuTime();
    tempo_reinsertion += (fimReinsertion - inicioreinsertion);
    return false;
  }
}

/*inline double calculateOrOpt2Cost(int i, int j, vector<int> &sequencia){
  double delta;

  if(i < j){
    delta = -matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i + 1]][sequencia[i + 2]] - matrizAdj[sequencia[j + 1]][sequencia[j + 2]] + matrizAdj[sequencia[i + 1]][sequencia[j + 2]] + matrizAdj[sequencia[i]][sequencia[j + 1]] + matrizAdj[sequencia[i - 1]][sequencia[i + 2]];
  }else{
    delta = -matrizAdj[sequencia[j - 1]][sequencia[j]] - matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i + 1]][sequencia[i + 2]] + matrizAdj[sequencia[j - 1]][sequencia[i]] + matrizAdj[sequencia[i + 1]][sequencia[j]] + matrizAdj[sequencia[i - 1]][sequencia[i + 2]];
  }
  return delta;
}*/

bool bestImprovementOrOpt2(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicioOropt2 = cpuTime();
  double melhorCusto = subseq_matrix[0][dimension].C;
  double custoTeste;
  //vector<int> copia = sequencia;
  int tamanho = dimension + 1;
  Subsequence sigma_1, sigma_2, sigma_3;

  int best_i = 0, best_j = 0;
  // cout << "tamanho sequencia: " << sequencia.size() << endl;

  for (int i = 1; i < tamanho - 2; i++){
    for (int j = 1; j < tamanho - 3; j++){
      // cout << "i: " << i << "j: " << j << endl;

      if (i != j)
      {
        /*for(int k = 0; k < copia.size() - 1; k++){
          cout << copia[k] << "->";
        }
          cout << copia.back() << endl;*/

        // copia.erase(copia.begin() + i, copia.begin() + i + 2);
        // copia.insert(copia.begin() + j, &sequencia[i], &sequencia[i] + 2);
        /*for(int k = 0; k < copia.size() - 1; k++){
          cout << copia[k] << "->";
        }
          cout << copia.back() << endl;*/
        if (i < j)
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+2][j]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[i][i+1]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[j+1][dimension-1]);    
        }
        else
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][j-1], subseq_matrix[i][i+1]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[j][i-1]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[i+2][dimension-1]);    
        }
        // delta = calculateOrOpt2Cost(i, j, sequencia);
        // custoTeste = custoDaSolucao(copia);

        /*if(custo + delta != custoTeste){
          cout << custo + delta << endl;
          cout << "custo teste: " << custoTeste << endl;
          cout << "i: " << i << "j: " << j << endl;
        }*/

        // cout << "Antes da copia" << endl;
        // copia = sequencia;

        if (sigma_3.C < melhorCusto)
        {
          melhorCusto = sigma_3.C;
          best_i = i;
          best_j = j;
        }
      }
    }
    // cout << "Teste apos o primeiro for" << endl;
  }
  // cout << "Teste 2" << endl;

  if (melhorCusto < subseq_matrix[0][dimension].C)
  {
    vector<int> reinserir(sequencia.begin() + best_i, sequencia.begin() + best_i + 2);
    sequencia.erase(sequencia.begin() + best_i, sequencia.begin() + best_i + 2);
    sequencia.insert(sequencia.begin() + best_j, reinserir.begin(), reinserir.end());
    UpdateAllSubseq(sequencia, subseq_matrix);
    // cout << "custo antes do Oropt2: " << custo << endl;
    //custo = custo + bestDelta;
    // cout << "custo apos o Oropt2: " << custo << endl;
    double fimOropt2 = cpuTime();
    tempo_orOpt2 += (fimOropt2 - inicioOropt2);
    return true;
  }
  else
  {
    double fimOropt2 = cpuTime();
    tempo_orOpt2 += (fimOropt2 - inicioOropt2);
    return false;
  }
}

/*inline double calculateOrOpt3Cost(int i, int j, vector<int> &sequencia){
  double delta;

  if(i < j){
    delta = -matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i + 2]][sequencia[i + 3]] - matrizAdj[sequencia[j + 2]][sequencia[j + 3]] + matrizAdj[sequencia[i - 1]][sequencia[i + 3]] + matrizAdj[sequencia[i]][sequencia[j + 2]] + matrizAdj[sequencia[i + 2]][sequencia[j + 3]];
  }else{
    delta = -matrizAdj[sequencia[j - 1]][sequencia[j]] - matrizAdj[sequencia[i - 1]][sequencia[i]] - matrizAdj[sequencia[i + 2]][sequencia[i + 3]] + matrizAdj[sequencia[j - 1]][sequencia[i]] + matrizAdj[sequencia[i + 2]][sequencia[j]] + matrizAdj[sequencia[i - 1]][sequencia[i + 3]];
  }
  return delta;
}*/

bool bestImprovementOrOpt3(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicioOropt3 = cpuTime();
  double melhorCusto = subseq_matrix[0][dimension].C;
  //vector<int> copia = sequencia;
  int tamanho = dimension + 1;
  Subsequence sigma_1, sigma_2, sigma_3;

  int best_i = 0, best_j = 0;
  // cout << "inicio" << endl;
  for (int i = 1; i < tamanho - 3; i++){
    for (int j = 1; j < tamanho - 4; j++){

      if (i != j)
      {
        /*copia.erase(copia.begin() + i, copia.begin() + i + 3);
        copia.insert(copia.begin() + j, &sequencia[i], &sequencia[i] + 3);*/
        if (i < j)
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][i-1], subseq_matrix[i+3][j]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[i][i+2]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[j+1][dimension-1]);      
        }
        else
        {
          sigma_1 = Subsequence::Concatenate(subseq_matrix[0][j-1], subseq_matrix[i][i+3]);
          sigma_2 = Subsequence::Concatenate(sigma_1, subseq_matrix[j][i-1]);
          sigma_3 = Subsequence::Concatenate(sigma_2, subseq_matrix[i+3][dimension-1]);     
        }
        // delta = calculateOrOpt3Cost(i, j, sequencia);

        // custoTeste = custoDaSolucao(copia);

        /*if(custo + delta != custoTeste){
          cout << custo + delta << endl;
          cout << custoTeste << endl;
          cout << "i: " << i << "j: " << j << endl;
        }*/

        // copia = sequencia;

        if (sigma_3.C < melhorCusto)
        {
          melhorCusto = sigma_3.C;
          best_i = i;
          best_j = j;
        }
      }
    }
  }

  if (melhorCusto < subseq_matrix[0][dimension].C)
  {
    vector<int> reinserir(sequencia.begin() + best_i, sequencia.begin() + best_i + 3);
    sequencia.erase(sequencia.begin() + best_i, sequencia.begin() + best_i + 3);
    sequencia.insert(sequencia.begin() + best_j, reinserir.begin(), reinserir.end());
    UpdateAllSubseq(sequencia, subseq_matrix);

    // cout << "custo antes do Oropt3: " << custo << endl;
    //custo = custo + bestDelta;
    // cout << "custo apos o Oropt3: " << custo << endl;

    double fimOropt3 = cpuTime();
    tempo_orOpt3 += (fimOropt3 - inicioOropt3);
    return true;
  }
  else
  {
    double fimOropt3 = cpuTime();
    tempo_orOpt3 += (fimOropt3 - inicioOropt3);
    return false;
  }
}

void buscaLocal(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  vector<int> NL = {1, 2, 3, 4, 5};
  bool improved = false;

  while (NL.empty() == false)
  {
    int n = rand() % NL.size();
    switch (NL[n])
    {
    case 1:
      improved = bestImprovementSwap(sequencia, subseq_matrix);
      break;
    case 2:
      improved = bestImprovement2Opt(sequencia, subseq_matrix);
      break;
    case 3:
      improved = bestImprovementReinsertion(sequencia, subseq_matrix);
      break;
    case 4:
      improved = bestImprovementOrOpt2(sequencia, subseq_matrix);
      break;
    case 5:
      improved = bestImprovementOrOpt3(sequencia, subseq_matrix);
      break;
    }

    if (improved)
    {
      NL = {1, 2, 3, 4, 5};
    }
    else
    {
      NL.erase(NL.begin() + n);
    }
  }
}

vector<int> perturbacao(vector<int> &sequencia, vector<vector<Subsequence>> &subseq_matrix)
{
  double inicioPert = cpuTime();
  int tamMaxSubSeq = ceil(((double)dimension) / 10);
  // cout << "tam max: " << tamMaxSubSeq << endl;
  int indice1, indice2;
  int tam1, tam2, indiceFinal1, indiceFinal2;
  int tamEntreSub; // tamanho entre subsequencias
  vector<int> copia = sequencia;

  indice1 = random(1, dimension - tamMaxSubSeq - tamMaxSubSeq - 1);
  // cout << "apos random 1" << endl;
  indiceFinal1 = random(indice1 + 1, indice1 + tamMaxSubSeq - 1); //+1 para que o indice final nao seja igual ao inicial
  // cout << "apos random 2" << endl;
  indice2 = random(indiceFinal1 + 1, dimension - tamMaxSubSeq);
  // cout << "apos random 3" << endl;
  indiceFinal2 = random(indice2 + 1, indice2 + tamMaxSubSeq - 1);
  // cout << "apos random 4" << endl;

  vector<int> subSeq1(copia.begin() + indice1, copia.begin() + indiceFinal1); // cria subsequencia 1 a partir do indice inicial e final 1
  vector<int> subSeq2(copia.begin() + indice2, copia.begin() + indiceFinal2); // cria subsequencia 2 a partir do indice inicial e final 2

  // cout << "subSeq1: " << endl;
  // printSolucao(subSeq1);
  // cout << "subSeq2: " << endl;
  // printSolucao(subSeq2);

  tamEntreSub = indice2 - indiceFinal1;
  tam1 = indiceFinal1 - indice1;
  tam2 = indiceFinal2 - indice2;

  // cout << "tam entre subsequencias: " << tamEntreSub << endl;
  // cout << "tam1: " << tam1 << endl;
  // cout << "tam2: " << tam2 << endl;

  copia.erase(copia.begin() + indice1, copia.begin() + indiceFinal1);    // apaga a subsequencia 1
  copia.insert(copia.begin() + indice1, subSeq2.begin(), subSeq2.end()); // insere a subsequencia 2 na posicao do indice 1
  // cout << "solucao antes do erase" << endl;
  // printSolucao(copia);
  copia.erase(copia.begin() + indice1 + tam2 + tamEntreSub, copia.begin() + indice1 + tam2 + tamEntreSub + tam2);
  // cout << "solucao depois do erase" << endl;
  // printSolucao(copia);
  copia.insert(copia.begin() + indice1 + tam2 + tamEntreSub, subSeq1.begin(), subSeq1.end());
  // cout << "solucao apos insert" << endl;
  // printSolucao(copia);

  // cout << "Solucao apos perturbacao: " << endl;
  // printSolucao(copia);

  double fimPert = cpuTime();

  tempo_perturbacao += (fimPert - inicioPert);

  return copia;
}

vector<int> ILS(int maxIter, int maxIterILS)
{
  vector<int> sequencia, best, bestOfAll;
  double custo = 0;
  double custoBest;
  double custoBestOfAll = INFINITY;

  for (int i = 0; i < maxIter; i++)
  {

    double alpha = (double)rand() / RAND_MAX;

    sequencia = construcao(alpha);

    vector<vector<Subsequence>> subseq_matrix(dimension + 1, vector<Subsequence>(dimension + 1));

    UpdateAllSubseq(sequencia, subseq_matrix);

    best = sequencia;

    // cout << "solucao inicial: " << endl;

    // printSolucao(sequencia);

    //custo = custoDaSolucao(sequencia);
    custoBest = subseq_matrix[0][dimension].C;

    // cout << "custo da solucao inicial: " << custo << endl;

    int iterILS = 0;

    while (iterILS <= maxIterILS)
    {
      buscaLocal(sequencia, subseq_matrix);
      if (subseq_matrix[0][dimension].C < custoBest)
      {
        // cout << "custo dentro do if: " << custo << endl;
        // cout << "dentro do if" << endl;
        custoBest = subseq_matrix[0][dimension].C;
        best = sequencia;
        // printSolucao(best);
        iterILS = 0;
      }
      // cout << "antes pert" << endl;
      sequencia = perturbacao(best, subseq_matrix);
      // cout << "sequencia apos pert" << endl;
      // printSolucao(sequencia);
      //custo = custoDaSolucao(sequencia);
      iterILS++;
    }
    if (custoBest < custoBestOfAll)
    {
      bestOfAll = best;
      custoBestOfAll = custoBest;
    }
  }
  return bestOfAll;
}

int main(int argc, char **argv)
{

  readData(argc, argv, &dimension, &matrizAdj);
  // printData();

  double somaTempos = 0;
  double somaValores = 0;

  int maxIter = 50;
  int maxIterILS;
  vector<int> solucaoIteracao;
  double valorIteracao;
  double valorTotal = numeric_limits<double>::infinity();

  if (dimension >= 150)
  {
    maxIterILS = dimension / 2.0;
  }
  else
  {
    maxIterILS = dimension;
  }

  for (int i = 0; i < 10; i++)
  {

    unsigned seed = time(0);
    // cout << seed << endl;
    srand(seed);

    double antes = cpuTime();

    solucaoIteracao = ILS(maxIter, maxIterILS);
    valorIteracao = custoDaSolucao(solucaoIteracao);

    double depois = cpuTime();

    double tempoTotal = depois - antes;

    somaTempos += tempoTotal;
    somaValores += valorIteracao;
  }

  if (valorTotal > valorIteracao)
  {
    valorTotal = valorIteracao;
    solucaoFinal = solucaoIteracao;
  }

  cout << "solucao best of all: " << endl;
  printSolucao(solucaoFinal);
  // custoBestOfAll = custoDaSolucao(bestOfAll);
  // cout << "Custo best of all: " << valorTotal << endl;

  cout << "Valor Medio da Solucao: " << (somaValores / 10) << endl;
  cout << "Tempo Medio de Execucao: " << (somaTempos / 10) << " (s)" << endl;

  printTime();

  // return bestOfAll;

  // perturbacao(sequencia);

  // printSolucao(sequencia);

  // bestImprovementSwap(sequencia, custo);

  // bestImprovement2Opt(sequencia, custo);

  // bestImprovementOrOpt2(sequencia, custo);

  // bestImprovementOrOpt3(sequencia, custo);

  // bestImprovementReinsertion(sequencia, custo);

  return 0;
}

void printData()
{
  cout << "dimension: " << dimension << endl;
  for (size_t i = 1; i <= dimension; i++)
  {
    for (size_t j = 1; j <= dimension; j++)
    {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}
