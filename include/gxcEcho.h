#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <array>

template <typename T>
void EchoVar(const std::string &VarName, const T &Var2Echo)
{ 
    std::cout << VarName << "=" << Var2Echo << std::endl;
    std::fstream logFile("gxmc.log", std::fstream::app);
    logFile << VarName << "=" << Var2Echo << std::endl;
    logFile.close();
}

void LogInitiate();
void EchoStr(const std::string &Str2Echo);
void EchoError(const std::string &Str2Echo);

void EchoMat(const std::string VarName, const std::vector<std::vector<double>> &Var2Echo);


template <typename T>
void EchoVarDebug(const T &Var2Echo){
    std::ofstream debugFile("debug.dat",std::ios_base::app);
    if (debugFile)
    {
        debugFile<<Var2Echo<<std::endl;
        debugFile.close();
    }
    else
        EchoStr("Unable to write debug file.");
    
}

void EchoVarDebugV3(const double Var2Echo[3]);

template <typename T>
void EchoArray(const std::string &VarName, int Num, T *Var2Echo)
{
    std::cout << VarName << "=";
    for (int i = 0; i < Num; i++)
    {
        std::cout << Var2Echo[i];
    }
    std::cout << std::endl;
}

void except(const bool &err, const std::string &errType);