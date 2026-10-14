#include "estimations.hpp"

using namespace std;
using namespace Eigen; 

double obtainMu0(){
    double mu0 = 0;
    for (int i = 1; i <= d; i++) {
        mu0+=log(v(0,i));
    }
    return mu0/d;
}

double obtainSigma0(){
    double sigma0 = 0;
    double mu0 = obtainMu0();

    for (int i = 1; i <= d; i++) {
        sigma0+=pow(log(v(0,i))-mu0,2);
    }
    return sigma0/d;
}

VectorXd arithmeticMean(){

    VectorXd mean(N+1);
    mean(0)=0;

    double sum;
    for(int j=1; j<=N; ++j){
        sum=0;
        for(int i=1; i<=d; ++i){
            sum+=x(i,j);
        }
        mean(j) = sum/d;
    }

    return mean;
}

VectorXd geometricMean(){
    
    VectorXd mean(N+1);
    mean(0)=0;

    for(int j=1; j<=N; ++j){
        double prod=1;
        for(int i=1; i<=d; ++i){
            prod*=x(i,j);
        }
        mean(j) = pow(prod,1.0/d);
    }

    return mean;
}

VectorXd vectorSigma2(){
        
    VectorXd sigma2(N+1);
    sigma2(0)=0;

    VectorXd ari(N+1);
    ari = arithmeticMean();
    VectorXd geo(N+1);
    geo = geometricMean();

    for(int j=1;j<=N;++j){
        sigma2(j) = 2*log(ari(j)/geo(j)); 
    }
    
    return sigma2;
}

VectorXd valuesPolymomialRegression(){
    
    VectorXd result(N);
    result(0)=0;

    VectorXd ari(N+1);
    ari = arithmeticMean();

    for(int j=1;j<=N-1;++j){
        double data = -log(ari(N)/ari(j) - 1); 
        if(!isnan(data)){
            result(j) = data;
        }
        else{
            result(j) = 0;
        }
    }

    return result;
}

double estimateSigma(){

    VectorXd sigma2(N+1);
    sigma2 = vectorSigma2();

    double sigma0 = obtainSigma0();

    VectorXd sigma_regression(N);

    for(int j=0;j<=N-1;++j){
        sigma_regression(j) =  sigma2(j+1) - sigma0; 
    }

    vector<double> regression = linearRegression(t_values,sigma_regression);

    function<double(double)> func = [regression](double x) { return regression[1] + regression[0] * x; };
    
    //graphPointsAndFunction(t_values, sigma_regression, func, 0, 50, "");    

    return regression[0];
}

VectorXd estimateEtaBeta(int degree){

    VectorXd y(N+1);
    y = valuesPolymomialRegression();

    VectorXd newY(N);
    newY = removeFirst(y);

    VectorXd newT(N);
    newT = removeLast(t_values);

    VectorXd poly(degree);
    poly = polynomialRegression(newT,newY,degree);

    function<double(double)> func = [poly](double x) {
        double result = 0.0;
        double power = 1.0;
        for (int i = 0; i < poly.size(); ++i) {
            result += poly[i] * power;
            power *= x;
        }
        return result;
    };

    //graphPointsAndFunction(newT,newY,func,0,50,"");

    return poly;
}
