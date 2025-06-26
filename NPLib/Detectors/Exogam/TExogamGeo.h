#ifndef TEXOGAMGEO_H
#define TEXOGAMGEO_H

#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

const double PI = 3.14159265358979323846;

class Vector3D {
public:

    Vector3D() : x(0.), y(0.), z(0.) {}
    Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    ~Vector3D(){}
    
    void setXYZ(double x_, double y_, double z_) {
        x = x_; y = y_; z = z_;
    }

    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    double scalarp(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3D vectorp(const Vector3D& other) const {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    void rotateAroundZ(double theta) {
        double x_ = x * cos(theta) - y * sin(theta);
        y = x * sin(theta) + y * cos(theta);
        x = x_;
    }
    void rotateAroundY(double phi) {
        double x_ = x * cos(phi) + z * sin(phi);
        z = -x * sin(phi) + z * cos(phi);
        x = x_;
    }
    
    void rotateAroundX(double theta) {
        double y_ = y * cos(theta) - z * sin(theta);
        z = y * sin(theta) + z * cos(theta);
        y = y_;
    }
    
    void rotateVector(double theta, double phi) {
        rotateAroundX(theta);
        rotateAroundY(phi);
    }
    
    double angle(const Vector3D& other) const {
        double ScalarProduct = this->scalarp(other);
        double NormsProduct = this->norm() * other.norm();
        return std::acos(ScalarProduct / NormsProduct);
    }
    
    void unit(){
        double norm = this->norm();
        x/=norm;
        y/=norm;
        z/=norm;
    }

    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    // Surcharge de l'opérateur de soustraction (-)
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }
    
    Vector3D operator*(double c) const{
        return Vector3D(c*x, c*y, c*z);
    }
    
    friend Vector3D operator*(double c, const Vector3D& v){
        return Vector3D(c*v.x, c*v.y, c*v.z);
    }
    
    
    Vector3D operator*(const float& c) const {
        return Vector3D(c*x, c*y, c*z);
    }

    void print(){
        std::cout << "Vector coordinates: (" << x << ", " << y << ", " << z << ")" << std::endl;
    }

    double X(){return x;}
    double Y(){return y;}
    double Z(){return z;}

    void clear(){x=0;y=0;z=0;}

    bool isnull(){return x==0 && y == 0 && z==0;}

    double theta() const {
       auto ZVector = Vector3D(0,0,1);
       return angle(ZVector);
    }

    double phi() const{
       auto XVector = Vector3D(1,0,0);
       auto ProjVector = Vector3D(x,y,0);
       if(y > 0)
        return ProjVector.angle(XVector);
       else
        return 2*PI - ProjVector.angle(XVector);
    }

    private:
        double x, y, z;
};

class TExogamGeo{

    protected:
        TExogamGeo();
    public:
        ~TExogamGeo(){};
        
        void SetBeam(const Vector3D& Beam_,const Vector3D& HitPosition_){Beam=Beam_; HitPosition=HitPosition_;};
        void SetRecoil(const Vector3D& Recoil_){Recoil=Recoil_;};
        void SetGammaInteractionPoint(const short& CloverNumber_,const short& CrystalNumber_,
        const short& OuterNumber_,const double& EGamma_); // Which outer + Energy gamma of the addbacked gamma-ray for better inter depth reconstruction
        // double GetGammaAngle(); // Returns angle between the recoil and the gamma interaction point
        // double GetDopplerE(); // Get Doppler Corrected E
        void ResetGamma(); // Resets Beam/HitPos/
        

        virtual Vector3D newton_raphson(double initial_guess,
        Vector3D& TargetPos, double tolerance, int max_iterations);
        double interaction_depth(const double Energy);
        virtual void setcrystal(Vector3D* SegmentPos) = 0; 
        virtual void setclover(Vector3D* SegmentPos) = 0; 

        void SetCloverNbr(const short& CloverNumber_);
        void SetCrystalNbr(const short& CrystalNumber_);
        void SetOuterNbr(const short& OuterNumberNumber_);
        void SetEGamma(const double& EGamma_);
        void SetInteractionDepth(const double& EGamma_);

    protected:
        double p_function(const double alpha, Vector3D& SegmentPos, const Vector3D& TargetPos, Vector3D& UnitVector);
        double equation(double alpha, Vector3D& SegmentPos, Vector3D& TargetPos, Vector3D& UnitVector);
        double equation_derivative(const double alpha, Vector3D& SegmentPos, Vector3D& TargetPos, Vector3D& UnitVector);
        Vector3D getouterpos(const double angle_alpha, const double InteractionDepth); 
        virtual double getdistance() = 0; 
        void ReadPhotonCS();
    
        
        Vector3D GammaDir;

        Vector3D Beam;
        Vector3D HitPosition;
        Vector3D Recoil;
        short CloverNumber;
        short CrystalNumber;
        short OuterNumber;
        double EGamma;
        double InteractionDepth;
        
        double GeDensity = 0.005323; // g/mm3
        std::map<double, double> Map_PhotonCS;
};
// Note: Each kind of ExogamGeometry is assumed to point at the center of the configuration
// If the target is offset from the center of the config, it is possible to precise a 
// TargetPos. If 2 ExogamClovers do not point on the same reference, it is recommended 
// to define 2 different geometries (but it is assumed that in a normal exp, the clovers 
// will all point to the same (0,0,0) point at the theoretical target position)
class TExogamStructure: public TExogamGeo{
    public:
        TExogamStructure(std::map<unsigned int, double> Distances_);
        ~TExogamStructure(){};
        
        
        
        void setclover(Vector3D* SegmentPos) override; 
        void setcrystal(Vector3D* SegmentPos) override; 
    private:
        std::map<unsigned int, double> Distances;
        double getdistance() override {return Distances[CloverNumber];}; 

};

class TExogamCartesian: public TExogamGeo{
    public:
        // clover nbr / crystal nbr / segment nbr / coordinates of the center of the segment
        TExogamCartesian(std::map<unsigned int, std::map<unsigned int,  std::map<unsigned int, Vector3D>>> Coordinates_);
        ~TExogamCartesian(){};
           
        void setclover(Vector3D* SegmentPos) override; 
        void setcrystal(Vector3D* SegmentPos) override; 
        Vector3D GetCloverCenter();
        Vector3D GetCrystalCenter(const unsigned int Crystal);

    private:
        std::map<unsigned int, std::map<unsigned int,  std::map<unsigned int, Vector3D>>> Coordinates;
        double getdistance() override; 
        double tolerance = 1; // tolerance is in mm
};
#endif 
