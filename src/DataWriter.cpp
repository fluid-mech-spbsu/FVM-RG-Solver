#include "datawriter.h"
#include <cstdio>
#include <cstdlib>
DataWriter::DataWriter(std::string pathName_)
{

    directory = pathName_;
    fs::path dataDir = directory / "data";
    fs::remove_all(dataDir);
    fs::create_directories(dataDir);
    fs::path infoFile = directory/"simInfo.txt";
    if(fileExist(infoFile.string()))
        std::remove(infoFile.string().c_str());
}

fs::path DataWriter::createTimeDirectory(double time)
{
    fs::path localDir = directory / "data" / std::to_string(time);
    fs::create_directories(localDir);
    return localDir;
}

void DataWriter::writeSimulationParam(string name, double value)
{
    ofstream info(directory/"simInfo.txt",std::ios::app);
    info<<name<<" : "<< value<< std::endl;
    info.close();
}

void DataWriter::writeSimulationParam(string name, string value)
{
    ofstream info(directory/"simInfo.txt",std::ios::app);
    info<<name<<" : "<< value<< std::endl;
    info.close();
}


void DataWriter::writeData(vector<macroParam> data, double time)
{
    fs::path localDir = createTimeDirectory(time);

    ofstream pressure(localDir/"pressure.txt",std::ios::out);
    ofstream velocity(localDir/"velocity.txt",std::ios::out);
    ofstream velocity_tau(localDir/"velocity_tau.txt",std::ios::out);
    ofstream velocity_normal(localDir/"velocity_normal.txt",std::ios::out);
    ofstream temp(localDir/"temp.txt",std::ios::out);
    ofstream density(localDir/"density.txt",std::ios::out);
//    ofstream densityO2(localDir/"densityO2.txt",std::ios::out);
//    ofstream densityO(localDir/"densityO.txt",std::ios::out);

//    ofstream e(localDir/"e.txt",std::ios::out);

    pressure<<"y"<<" "<<"p"<<endl;
    velocity<<"y"<<" "<<"v"<<endl;
    velocity_tau<<"y"<<" "<<"v_t"<<endl;
    velocity_normal<<"y"<<" "<<"v_n"<<endl;
    temp<<"y"<<" "<<"T"<<endl;
    density<<"y"<<" "<<"rho"<<endl;
//    densityO2<<"y"<<" "<<"rho_O2"<<endl;
//    densityO<<"y"<<" "<<"rho_O"<<endl;

    // pressure<< dh/2 <<" "<<data[0].pressure<<endl;
    // velocity<< dh/2 <<" "<<(data[0].velocity + data[1].velocity) / 2. <<endl;
    // velocity_tau<< dh/2 <<" "<<(data[0].velocity_tau + data[1].velocity_tau) / 2. <<endl;
    // velocity_normal<< dh/2 <<" "<<(data[0].velocity_normal + data[1].velocity_normal) / 2. <<endl;
    // temp<< dh/2 <<" "<< (data[0].temp + data[1].temp) / 2. <<endl;
    // density<< dh/2 <<" " <<data[0].density <<endl;  

    pressure<< 0 <<" "<<data[0].pressure<<endl;
    velocity<< 0 <<" "<<(data[0].velocity + data[1].velocity) / 2. <<endl;
    velocity_tau<< 0 <<" "<<(data[0].velocity_tau + data[1].velocity_tau) / 2. <<endl;
    velocity_normal<< 0 <<" "<<(data[0].velocity_normal + data[1].velocity_normal) / 2. <<endl;
    temp<< 0 <<" "<< (data[0].temp + data[1].temp) / 2. <<endl;
    density<< 0 <<" " <<data[0].density <<endl;

    for(size_t i = 1; i < data.size() - 1; i++)
    {
        pressure<<dh*i - dh/2<<" "<<data[i].pressure<<endl;
        // pressure<<dh*i<<" "<<data[i].pressure<<endl;
        velocity<<dh*i - dh/2<<" "<<data[i].velocity<<endl;
        velocity_tau<<dh*i - dh/2<<" "<<data[i].velocity_tau<<endl;
        velocity_normal<<dh*i - dh/2<<" "<<data[i].velocity_normal<<endl;
        temp<<dh*i - dh/2<<" "<<data[i].temp<<endl;
        density<<dh*i - dh/2<<" "<<data[i].density<<endl;
//        if(data[i].mixture.NumberOfComponents == 2)
//        {
//            densityO2<<dh*i<<" "<<data[i].densityArray[0]<<endl;
//            densityO<<dh*i<<" "<<data[i].densityArray[1]<<endl;
//        }

    }

    // double dh_last = dh * (data.size()-1) - dh/2;
    double dh_last = dh * (data.size()-1) - dh;
    size_t s = data.size()-1;
    pressure<<dh_last<<" "<<data[s].pressure<<endl;
    velocity<<dh_last<<" "<<(data[s].velocity + data[s-1].velocity) / 2.<<endl;
    velocity_tau<<dh_last<<" "<<(data[s].velocity_tau + data[s-1].velocity_tau) / 2.<<endl;
    velocity_normal<<dh_last<<" "<<(data[s].velocity_normal + data[s-1].velocity_normal) / 2.<<endl;
    temp<<dh_last<<" "<<(data[s].temp + data[s-1].temp) / 2.<<endl;
    density<<dh_last<<" "<<data[s].density<<endl;

    pressure.close();
    velocity.close();
    temp.close();
    velocity_normal.close();
    velocity_tau.close();
}

void DataWriter::writeData(vector<macroParam> data, macroParam dataDown, macroParam dataUp, double time)
{
    fs::path localDir = createTimeDirectory(time);

    ofstream pressure(localDir/"pressure.txt",std::ios::out);
    ofstream velocity(localDir/"velocity.txt",std::ios::out);
    ofstream velocity_tau(localDir/"velocity_tau.txt",std::ios::out);
    ofstream velocity_normal(localDir/"velocity_normal.txt",std::ios::out);
    ofstream temp(localDir/"temp.txt",std::ios::out);
    ofstream density(localDir/"density.txt",std::ios::out);

    pressure<<"y"<<" "<<"p"<<endl;
    velocity<<"y"<<" "<<"v"<<endl;
    velocity_tau<<"y"<<" "<<"v_t"<<endl;
    velocity_normal<<"y"<<" "<<"v_n"<<endl;
    temp<<"y"<<" "<<"T"<<endl;
    density<<"y"<<" "<<"rho"<<endl;

    pressure<<0<<" "<<dataDown.pressure<<endl;
    velocity<<0<<" "<<dataDown.velocity<<endl;
    velocity_tau<<0<<" "<<dataDown.velocity_tau<<endl;
    velocity_normal<<0<<" "<<dataDown.velocity_normal<<endl;
    temp<<0<<" "<<dataDown.temp<<endl;
    density<<0<<" "<<dataDown.density<<endl;

    for(size_t i = 1; i < data.size()-1; i++)
    {
        pressure<<dh*i<<" "<<data[i].pressure<<endl;
        velocity<<dh*i<<" "<<data[i].velocity<<endl;
        velocity_tau<<dh*i<<" "<<data[i].velocity_tau<<endl;
        velocity_normal<<dh*i<<" "<<data[i].velocity_normal<<endl;
        temp<<dh*i<<" "<<data[i].temp<<endl;
        density<<dh*i<<" "<<data[i].density<<endl;
    }
    double dh_last = dh * (data.size()-1);
    pressure<<dh_last<<" "<<dataUp.pressure<<endl;
    velocity<<dh_last<<" "<<dataUp.velocity<<endl;
    velocity_tau<<dh_last<<" "<<dataUp.velocity_tau<<endl;
    velocity_normal<<dh_last<<" "<<dataUp.velocity_normal<<endl;
    temp<<dh_last<<" "<<dataUp.temp<<endl;
    density<<dh_last<<" "<<dataUp.density<<endl;

    pressure.close();
    velocity.close();
    temp.close();
    velocity_normal.close();
    velocity_tau.close();
}

void DataWriter::setDelta_h(double dh_)
{
    dh = dh_;
}

bool DataWriter::fileExist(const string &name)
{
    ifstream f(name.c_str());
    return f.good();
}


bool DataReader::read()
{
    if(!fillDataVector(pres,"/pressure.txt"))
        return 0;
    if(!fillDataVector(vel,"/velocity.txt"))
        return 0;
    if(!fillDataVector(temp,"/temp.txt"))
        return 0;
    if(!fillDataVector(density,"/density.txt"))
        return 0;
    if(!fillDataVector(vel_tau,"/velocity_tau.txt"))
        return 0;
    if(!fillDataVector(vel_normal,"/velocity_normal.txt"))
        return 0;
    return 1;
}

void DataReader::getPoints(vector<macroParam> &points)
{
    points.clear();
    for(size_t i = 0; i < pres.size();i++)
    {
        macroParam tmp;
        tmp.densityArray = {density[i]};
        tmp.fractionArray = {1};
        tmp.density      = density[i];
        tmp.pressure     = pres[i];
        tmp.velocity_tau     =  vel_tau[i];
        tmp.velocity_normal  = vel_normal[i];
        tmp.velocity = vel[i];
        tmp.temp         = temp[i];
        tmp.gas         ="Ar";
        points.push_back(tmp);
    }
    return;
}

bool DataReader::fillDataVector(vector<double> &data, string dataFileName)
{
    std::ifstream file(pathName + dataFileName, std::ios::in); // окрываем файл для чтения
    if (!file.is_open())
    {
        cout<<"WARNING: no" << dataFileName <<" file to read"<<endl;
        return 0;
    }
    double h, value;
    string y,var;
    file>>y>>var;
    while (file >> h >> value)
    {
        data.push_back(value);
    }
    file.close();     // закрываем файл
    return 1;
}
