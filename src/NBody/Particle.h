/*! \file Particle.h
 *  \brief header file for the \ref NBody::Particle and derived classes

*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <NBodyMath.h>
#include <SwiftParticle.h>
#include <map>
#include <vector>
#include <memory>

#ifdef USEBOOSTMPI
#include <boost/mpi.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#endif

using namespace std;
using namespace Math;
#ifdef SWIFTINTERFACE
using namespace Swift;
#endif
namespace NBody
{


/// \name Particle Memory Flags
/// \brief These are define flags that alter the memory allocation of particles
//@{
///use long ints to store particle ids
///if NOMASS is defined, don't actually store a mass
//#define NOMASS
///this is the mass value that is returned if no value is actually stored
#define MASSVAL 1.0
///use floats to store positions and velocities regardless of the precision used elsewhere
#ifdef LOWPRECISIONPOS
typedef float DoublePos_t;
#else
typedef Double_t DoublePos_t;
#endif
///define ID (index) type, whether signed or unsigned Int_t
#ifdef PARTICLEUIDS
typedef UInt_t PARTIDTYPE;
#else
typedef Int_t PARTIDTYPE;
#endif

///define permenant (particle) ID type, whether signed or unsigned Int_t
#ifdef PARTICLEUPIDS
typedef UInt_t PARTPIDTYPE;
#else
typedef Int_t PARTPIDTYPE;
#endif
//@}


    /*!
    \class NBody::HydroProperties
    \brief A simple class to store hydrodynamic quantities
    */
    class HydroProperties
    {
        protected:
            ///store properties like SphDen, SphPressure, Volume, SelfEnergy, etc
            map<string, float> InternalProperties;
            map<string, float> Chemistry;
            map<string, float> Feedback;
            map<string, float> ChemistryProduction;
        public:
            HydroProperties(){};
            HydroProperties(const HydroProperties&) = default;
            HydroProperties(HydroProperties&&) = default;
            HydroProperties& operator=(const HydroProperties&) = default;
            HydroProperties& operator=(HydroProperties&&) = default;
            bool operator==(const HydroProperties &h) const
            {
                int ival = 1;
                ival *= (InternalProperties == h.InternalProperties);
                ival *= (Chemistry == h.Chemistry);
                ival *= (Feedback == h.Feedback);
                ival *= (ChemistryProduction == h.ChemistryProduction);
                return bool(ival);
            };
            ~HydroProperties(){};

            float GetInternalProperties(string &f) {return InternalProperties[f];}
            float GetChemistry(string &f) {return Chemistry[f];}
            float GetFeedback(string &f) {return Feedback[f];}
            float GetChemistryProduction(string &f) {return ChemistryProduction[f];}

            map<string, float> GetInternalProperties() {return InternalProperties;}
            map<string, float> GetChemistry() {return Chemistry;}
            map<string, float> GetFeedback() {return Feedback;}
            map<string, float> GetChemistryProduction() {return ChemistryProduction;}

            void SetInternalProperties(string &f, float value) {InternalProperties[f]=value;}
            void SetChemistry(string &f, float value) {Chemistry[f]=value;}
            void SetFeedback(string &f, float value) {Feedback[f]=value;}
            void SetChemistryProduction(string &f, float value) {ChemistryProduction[f]=value;}

            void SetInternalProperties(map<string, float> m) {InternalProperties=m;}
            void SetChemistry(map<string, float> m) {Chemistry=m;}
            void SetFeedback(map<string, float> m) {Feedback=m;}
            void SetChemistryProduction(map<string, float> m) {ChemistryProduction=m;}
            //example eagle chemistry is {"Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen","Neon", "Magnesium", "Silicon", "Iron"};
    };

    /*!
    \class NBody::StarProperties
    \brief A simple class to store stellar quantities
    */
    class StarProperties
    {
        protected:
            map<string, float> InternalProperties;
            map<string, float> Chemistry;
            map<string, float> Feedback;
            map<string, float> ChemistryProduction;
        public:
            StarProperties(){};
            StarProperties(const StarProperties&) = default;
            StarProperties(StarProperties&&) = default;
            StarProperties& operator=(const StarProperties&) = default;
            StarProperties& operator=(StarProperties&&) = default;
            bool operator==(const StarProperties &s) const
            {
                int ival = 1;
                ival *= (InternalProperties == s.InternalProperties);
                ival *= (Chemistry == s.Chemistry);
                ival *= (Feedback == s.Feedback);
                ival *= (ChemistryProduction == s.ChemistryProduction);
                return bool(ival);
            };
            ~StarProperties(){};

            float GetInternalProperties(string &f) {return InternalProperties[f];}
            float GetChemistry(string &chem) {return Chemistry[chem];};
            float GetFeedback(string &f) {return Feedback[f];};
            float GetChemistryProduction(string &f) {return ChemistryProduction[f];};

            map<string, float> GetInternalProperties() {return InternalProperties;}
            map<string, float> GetChemistry() {return Chemistry;}
            map<string, float> GetFeedback() {return Feedback;}
            map<string, float> GetChemistryProduction() {return ChemistryProduction;}

            void SetInternalProperties(string &f, float value) {InternalProperties[f]=value;}
            void SetChemistry(string &chem, float value) {Chemistry[chem]=value;};
            void SetFeedback(string &f, float value) {Feedback[f]=value;};
            void SetChemistryProduction(string &f, float value) {ChemistryProduction[f]=value;};

            void SetInternalProperties(map<string, float> m) {InternalProperties=m;}
            void SetChemistry(map<string, float> m) {Chemistry=m;}
            void SetFeedback(map<string, float> m) {Feedback=m;}
            void SetChemistryProduction(map<string, float> m) {ChemistryProduction=m;}
    };
    /*!
    \class NBody::BHProperties
    \brief A simple class to store black hole  quantities
    */
    class BHProperties
    {
        protected:
            map<string, float> InternalProperties;
            map<string, float> Chemistry;
            map<string, float> Feedback;
            map<string, float> ChemistryProduction;
            map<string, float> AccretedMassChannel;
        public:
            BHProperties(){};
            BHProperties(const BHProperties&) = default;
            BHProperties(BHProperties&&) = default;
            BHProperties& operator=(const BHProperties&) = default;
            BHProperties& operator=(BHProperties&&) = default;
            bool operator==(const BHProperties &b) const
            {
                int ival = 1;
                ival *= (InternalProperties == b.InternalProperties);
                ival *= (Chemistry == b.Chemistry);
                ival *= (Feedback == b.Feedback);
                ival *= (ChemistryProduction == b.ChemistryProduction);
                ival *= (AccretedMassChannel == b.AccretedMassChannel);
                return bool(ival);
            };
            ~BHProperties(){};

            float GetInternalProperties(string &f) {return InternalProperties[f];}
            float GetChemistry(string &chem) {return Chemistry[chem];};
            float GetFeedback(string &f) {return Feedback[f];};
            float GetChemistryProduction(string &f) {return ChemistryProduction[f];};
            float GetAccretedMassChannel(string &f) {return AccretedMassChannel[f];};

            map<string, float> GetInternalProperties() {return InternalProperties;}
            map<string, float> GetChemistry() {return Chemistry;}
            map<string, float> GetFeedback() {return Feedback;}
            map<string, float> GetChemistryProduction() {return ChemistryProduction;}
            map<string, float> GetAccretedMassChannel() {return AccretedMassChannel;}

            void SetInternalProperties(string &f, float value) {InternalProperties[f]=value;}
            void SetChemistry(string &f, float value) {Chemistry[f]=value;};
            void SetFeedback(string &f, float value) {Feedback[f]=value;};
            void SetChemistryProduction(string &f, float value) {ChemistryProduction[f]=value;};
            void SetAccretedMassChannel(string &f, float value) {AccretedMassChannel[f]=value;};

            void SetInternalProperties(map<string, float> m) {InternalProperties=m;}
            void SetChemistry(map<string, float> m) {Chemistry=m;}
            void SetFeedback(map<string, float> m) {Feedback=m;}
            void SetChemistryProduction(map<string, float> m) {ChemistryProduction=m;}
            void SetAccretedMassChannel(map<string, float> m) {AccretedMassChannel=m;}
    };

    /*!
    \class NBody::ExtraDMProperties
    \brief A simple class to store extra dark matter (gravitaional particle) properties
    */
    class ExtraDMProperties
    {
        protected:
            map<string, float> ExtraProperties;
        public:
            ExtraDMProperties(){};
            ExtraDMProperties(const ExtraDMProperties&) = default;
            ExtraDMProperties(ExtraDMProperties&&) = default;
            ExtraDMProperties& operator=(const ExtraDMProperties&) = default;
            ExtraDMProperties& operator=(ExtraDMProperties&&) = default;
            bool operator==(const ExtraDMProperties &d) const
            {
                int ival = 1;
                ival *= (ExtraProperties == d.ExtraProperties);
                return bool(ival);
            };
            ~ExtraDMProperties() = default;

            float GetExtraProperties(string &f){return ExtraProperties[f];}
            map<string, float> GetExtraProperties() {return ExtraProperties;}
            void SetExtraProperties(string &f, float value) {ExtraProperties[f]=value;};
            void SetExtraProperties(map<string, float> m) {ExtraProperties=m;}
    };
/*!
    \class NBody::Particle
    \brief A simple n-body particle class.

    The class is for a simple n-body particle with mass, position, velocity, id, type and density information.
    Note that here particles can be compiled to store mass or assume all particles have the same mass given by a
    MASSVAL definition in \ref Particle.h by defining the flag NOMASS. Particles can also be compiled to have long ids
    stored using an unsigned int or unsigned long it.
*/
    class Particle
    {
        protected:
        ///mass
#ifndef NOMASS
        Double_t mass;
#endif
#ifdef LOWPRECISIONPOS
        /// x, y, z
        DoublePos_t position[3];
        /// vx, vy, vz
        DoublePos_t velocity[3];
#else
        /// x, y, z
        Double_t position[3];
        /// vx, vy, vz
        Double_t velocity[3];
#endif
        ///\name identifiers
        //@{
        //extra particle identifier, useful as it is not changed by KDTree
        PARTPIDTYPE pid;
        PARTIDTYPE id;
        //int gid;
        int type;
        //@}
        ///density
        Double_t rho;
        ///potential
        Double_t phi;
#ifdef SWIFTINTERFACE
        Double_t gravityphi;
        int swifttask;
        int swiftindex;
#endif
        ///For hydrodynamical quantities
        //@{
        ///if gas flag is set, then particle also can have sph based quantities. if star flag is set,
        ///then particle also can have star properties such as age and metallicity
#ifdef GASON
        ///self energy
        DoublePos_t u;
        ///sph based density
        DoublePos_t sphden;
        unique_ptr<HydroProperties> hydro;
#endif
#ifdef STARON
        ///stellar age
        DoublePos_t tage;
        unique_ptr<StarProperties> star;
#endif
#if defined (GASON) && (STARON)
        ///metallicity
        DoublePos_t zmet;
        ///star formation rate of gas
        DoublePos_t sfr;
#endif
#if (defined(GASON) && defined(GASEXTRA)) || (defined(GASON) && defined(SWIFTINTERFACE))
        ///store entropy, useful for searching for shocks and related entropy structures
        DoublePos_t entropy;
        ///store temperature
        DoublePos_t temperature;
#endif
        //@}

#ifdef EXTRAINPUTINFO
        ///Variables to store extra info related to where particle read from some input data
        ///used by VELOCIraptor group finding
        //@{
        ///file id of input file containing particle
        Int_t inputFileID;
        ///index of particle in input file
        Int_t inputIndexInFile;
        //@}
#endif
#ifdef EXTRAFOFINFO
        ///Variables to store VELOCIraptor group finding info
        //@{
        Int_t GroupID;
        Int_t ParentGroupID;
        Int_t SOGroupID;
        //@}
#endif

#ifdef BHON
        unique_ptr<BHProperties> bh;
#endif

#ifdef EXTRADMON
        unique_ptr<ExtraDMProperties> dm;
#endif

        public:

        /// \name Constructors & Destructors
        //@{
        Particle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, PARTPIDTYPE PID=0);
        Particle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, PARTPIDTYPE PID=0);
        Particle(const Particle &p);
        Particle(Particle &&p) = default;
#ifdef SWIFTINTERFACE
        Particle(const struct gpart &p,  double lscale, double vscale, double mscale, double uscale, bool icosmological=true, double scalefactor=1.0, double littleh=1.0);
        Particle(const struct swift_vel_part &p);
#endif
        Particle(std::istream &F);
        //No dynamic allocation, thus destructor not needed.
        ~Particle(){};
        //@}


        /// \name Overloaded operators
        //@{
        Particle& operator=(const Particle &part);
        Particle& operator=(Particle&&) = default;
        bool operator==(const Particle &p) const
        {
            int ival=1;
            ival*=((position[0]==p.position[0])&&(position[1]==p.position[1])&&(position[2]==p.position[2])&&
                (velocity[0]==p.velocity[0])&&(velocity[1]==p.velocity[1])&&(velocity[2]==p.velocity[2])&&
                (id==p.id)&&(type==p.type)&&(rho==p.rho)&&(phi==p.phi)&&
                (pid==p.pid));
#ifdef SWIFTINTERFACE
            ival*=((gravityphi==p.gravityphi)&&(swifttask==p.swifttask)&&(swiftindex==p.swiftindex));
#endif
#ifndef NOMASS
            ival*=((mass==p.mass));
#endif
#ifdef GASON
            ival*=((u==p.u)&&(sphden==p.sphden));
#endif
#ifdef STARON
            ival*=((tage==p.tage));
#endif
#if defined(GASON) && defined(STARON)
            ival*=((zmet==p.zmet)&&(sfr==p.sfr));
#endif

#if (defined(GASON) && defined(GASEXTRA)) || (defined(GASON) && defined(SWIFTINTERFACE))
            ival*=((entropy==p.entropy));
            ival*=((temperature==p.temperature));
#endif
#ifdef EXTRAINPUTINFO
            ival*=((inputFileID==p.inputFileID));
            ival*=((inputIndexInFile==p.inputIndexInFile));
#endif
#ifdef EXTRAFOFINFO
            ival*=((GroupID==p.GroupID));
            ival*=((ParentGroupID==p.ParentGroupID));
            ival*=((SOGroupID==p.SOGroupID));
#endif
#ifdef GASON
            if (hydro && p.hydro) ival *= (*hydro == *p.hydro);
            else ival *=0;
#endif
#ifdef STARON
            if (star && p.star) ival *= (*star == *p.star);
            else ival *=0;
#endif
#ifdef BHON
            if (bh && p.bh) ival *= (*bh == *p.bh);
            else ival *=0;
#endif
#ifdef EXTRADMON
            if (dm && p.dm) ival *= (*dm == *p.dm);
            else ival *=0;
#endif
            return ival;
        }
        bool operator!=(const Particle &p) const {return !((*this)==p);}
        // Inserter: prints out the particle in ascii format,
        // e.g., cout << part;
        friend ostream &operator << (ostream &outs, const Particle &p);
        // Extractor: sets particle attributes from a stream,
        // e.g., cin >> part;
        //friend istream &operator >> (istream &ins, Particle &p);
        //@}

        /// \name IO methods
        //@{
        //// Read a particle from an stream in binary format
        //void Read(std::istream &F);
        /// Write a particle to a stream in binary format
        void Write(std::ostream &F) const;
        //@}

        /// \name Get & Set methods
        //@{
#ifdef NOMASS
        Double_t GetMass() const { return MASSVAL; }
        void SetMass(const Double_t &Mass) {}
#else
        Double_t GetMass() const { return mass; }
        void SetMass(const Double_t &Mass) {mass=Mass;}
#endif

#ifndef LOWPRECISIONPOS
        const Double_t* GetPosition() const { return position; }
#else
        const DoublePos_t* GetPosition() const { return position; }
#endif
        Double_t GetPosition(const int &i) const { return position[i]; }
        Double_t X() const { return position[0]; }
        Double_t Y() const { return position[1]; }
        Double_t Z() const { return position[2]; }
        void SetPosition(const int &i, const Double_t &x){ position[i] = x; }
        void SetPosition(const Double_t* p)
        {
        position[0] = p[0];
        position[1] = p[1];
        position[2] = p[2];
        }
        void SetPosition(const Double_t &x, const Double_t &y, const Double_t &z)
        {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        }

#ifndef LOWPRECISIONPOS
        const Double_t *GetVelocity() const { return velocity; }
#else
        const DoublePos_t* GetVelocity() const { return velocity; }
#endif
        Double_t GetVelocity(const int &i) const { return velocity[i]; }
        Double_t Vx() const { return velocity[0]; }
        Double_t Vy() const { return velocity[1]; }
        Double_t Vz() const { return velocity[2]; }
        void SetVelocity(const int &i, const Double_t &x){ velocity[i] = x; }
        void SetVelocity(const Double_t* v)
        {
        velocity[0] = v[0];
        velocity[1] = v[1];
        velocity[2] = v[2];
        }
        void SetVelocity(const Double_t &vx, const Double_t &vy, const Double_t &vz)
        {
        velocity[0] = vx;
        velocity[1] = vy;
        velocity[2] = vz;
        }
        Double_t GetPhase(const int &i) const {
            if (i<3) return position[i];
            else return velocity[i-3];
        }

        void SetPhase(const int &i, const Double_t &x){
            if (i<3) position[i]=x;
            else velocity[i-3]=x;
        }

        PARTPIDTYPE GetPID() const {return pid;}
        void SetPID(const PARTPIDTYPE &i){pid=i;}
        PARTIDTYPE GetID() const {return id;}
        void SetID(const PARTIDTYPE &i){id=i;}
        int GetType() const {return type;}
        void SetType(int i){type=i;}
        Double_t GetDensity() const {return rho;}
        void SetDensity(const Double_t &Rho){rho=Rho;}
        Double_t GetPotential() const {return phi;}
        void SetPotential(const Double_t &Phi){phi=Phi;}
#ifdef SWIFTINTERFACE
        Double_t GetGravityPotential() {return gravityphi;}
        void SetGravityPotential(const Double_t &Phi){gravityphi=Phi;}
        int GetSwiftIndex() {return swiftindex;}
        void SetSwiftIndex(const int &index){swiftindex=index;}
        int GetSwiftTask() {return swifttask;}
        void SetSwiftTask(const int &task){swifttask=task;}
#endif
#ifdef GASON
        Double_t GetU() const {return u;}
        void SetU(const Double_t &U){u=U;}
        Double_t GetSPHDen() const {return sphden;}
        void SetSPHDen(const Double_t &SPHden){sphden=SPHden;}
        ///\todo having temperature and entropy functions are not trivial
        ///because need to include eos of gas, metallicity, units, etc.
#endif
#ifdef STARON
        Double_t GetTage() const {return tage;}
        void SetTage(const Double_t &Tage){tage=Tage;}
#endif
#if defined (GASON) && (STARON)
        Double_t GetZmet() const {return zmet;}
        void SetZmet(const Double_t &Zmet){zmet=Zmet;}
        Double_t GetSFR() const {return sfr;}
        void SetSFR(const Double_t &SFR){sfr=SFR;}
#endif

#if (defined(GASON) && defined(GASEXTRA)) || (defined(GASON) && defined(SWIFTINTERFACE))
        Double_t GetEntropy() const {return entropy;}
        void SetEntropy(const Double_t &Entropy) {entropy=Entropy;}
        Double_t GetTemperature() const {return temperature;}
        void SetTemperature(const Double_t &Temperature) {temperature=Temperature;}
#endif

#ifdef EXTRAINPUTINFO
        ///Sets and Gets for ExtendedOutput variables
        void SetInputFileID(const Int_t &i) {inputFileID = i;}
        Int_t GetInputFileID() const {return inputFileID;}
        void SetInputIndexInFile(const Int_t &i) {inputIndexInFile = i;}
        Int_t GetInputIndexInFile() const {return inputIndexInFile;}
#endif
#ifdef EXTRAFOFINFO
        ///Sets and Gets for ExtendedOutput variables
        void SetGroupID(const Int_t &i) {GroupID = i;}
        Int_t GetGroupID() const {return GroupID;}
        void SetParentGroupID(const Int_t &i) {ParentGroupID = i;}
        Int_t GetParentGroupID() const {return ParentGroupID;}
        void SetSOGroupID(const Int_t &i) {SOGroupID = i;}
        Int_t GetSOGroupID() const {return SOGroupID;}
#endif

#ifdef GASON
        bool HasHydroProperties() { return bool(hydro); };
        void InitHydroProperties() { hydro.reset(new HydroProperties()); };
        HydroProperties& GetHydroProperties() {return *hydro;};
        void SetHydroProperties(const HydroProperties &value) {
            hydro.reset(new HydroProperties(value));
        };
        void SetHydroProperties() {hydro.reset(nullptr);};
        ///function call that releases ownership of the pointer
        HydroProperties* ReleaseHydroProperties() { return hydro.release(); };
        ///function that releases a pointer and sets that pointer to NULL.
        ///To be used in specific circumstances such as MPI byte copies
        ///as will otherwise lead to memory leaks. NOT IDEAL
        void NullHydroProperties()
        {
            HydroProperties *h = hydro.release();
            h=nullptr;
            hydro.reset(nullptr);
        };
#endif
#ifdef STARON
        bool HasStarProperties() { return bool(star); };
        void InitStarProperties() { star.reset(new StarProperties()); };
        StarProperties& GetStarProperties() {return *star;};
        void SetStarProperties(const StarProperties &value) {
            star.reset(new StarProperties(value));
        };
        void SetStarProperties() {star.reset(nullptr);};
        ///function call that releases ownership of the pointer
        StarProperties* ReleaseStarProperties() { return star.release(); };
        ///function that releases a pointer and sets that pointer to NULL.
        ///To be used in specific circumstances such as MPI byte copies
        ///as will otherwise lead to memory leaks. NOT IDEAL
        void NullStarProperties()
        {
            StarProperties *s = star.release();
            s=nullptr;
            star.reset(nullptr);
        };
#endif
#ifdef BHON
        bool HasBHProperties() { return bool(bh); };
        void InitBHProperties() { bh.reset(new BHProperties()); };
        BHProperties& GetBHProperties() {return *bh;};
        void SetBHProperties(const BHProperties &value) {
            bh.reset(new BHProperties(value));
        };
        void SetBHProperties() {bh.reset(nullptr);};
        ///function call that releases ownership of the pointer
        BHProperties* ReleaseBHProperties() { return bh.release(); };
        ///function that releases a pointer and sets that pointer to NULL.
        ///To be used in specific circumstances such as MPI byte copies
        ///as will otherwise lead to memory leaks. NOT IDEAL
        void NullBHProperties()
        {
            BHProperties *b = bh.release();
            b=nullptr;
            bh.reset(nullptr);
        };
#endif

#ifdef EXTRADMON
        bool HasExtraDMProperties() { return bool(dm); };
        void InitExtraDMProperties() { dm.reset(new ExtraDMProperties()); };
        ExtraDMProperties& GetExtraDMProperties() {return *dm;};
        void SetExtraDMProperties(const ExtraDMProperties &value) {
            dm.reset(new ExtraDMProperties(value));
        };
        void SetExtraDMProperties() {dm.reset(nullptr);};
        ///function call that releases ownership of the pointer
        ExtraDMProperties* ReleaseExtraDMProperties() { return dm.release(); };
        ///function that releases a pointer and sets that pointer to NULL.
        ///To be used in specific circumstances such as MPI byte copies
        ///as will otherwise lead to memory leaks. NOT IDEAL
        void NullExtraDMProperties()
        {
            ExtraDMProperties *d = dm.release();
            d=nullptr;
            dm.reset(nullptr);
        };
#endif
        //@}

        /// \name Other useful functions
        //@{
        ///scale positions by a quantity
        void ScalePositions(Double_t &x){position[0]*=x;position[1]*=x;position[2]*=x;}
        ///scale positions by a quantity
        void ScaleVelocity(Double_t &x){velocity[0]*=x;velocity[1]*=x;velocity[2]*=x;}
        ///scale positions/velocity by a quantity
        void ScalePhase(Double_t &x, Double_t &v){position[0]*=x;position[1]*=x;position[2]*=x;velocity[0]*=v;velocity[1]*=v;velocity[2]*=v;}
        /// Radius of particle
        Double_t Radius() const {return sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);}
        ///Radius squared
        Double_t Radius2() const{return position[0] * position[0]+position[1]*position[1]+position[2]*position[2];}
        /// Spherical coordinate theta (angle between the position vector and the z-axis)
        Double_t Theta() const
        {
            Double_t r=Radius();
            if (r>0) return acos(position[2] / Radius());
            else return 0;
        }
        /// Spherical coordinate phi (azimuthal angle)
        Double_t Phi() const
        {
            if (position[1]>=0&&position[0]>=0)
            return atan2(position[1], position[0]);
            else if (position[1]>0&&position[0]<0)
            return -atan2(position[1], -position[0])+1.0*M_PI;
            else if (position[1]<0&&position[0]<0)
            return atan2(-position[1], -position[0])+1.0*M_PI;
            else
            return -atan2(-position[1], position[0])+2.0*M_PI;
        }
        /// Cylindrical radius (with symmetry axis along z)
        Double_t CylRadius() const
        {
            return sqrt(position[0] * position[0] + position[1] * position[1]);
        }
        /// Absolute velocity of particle
        Double_t AbsoluteVelocity() const
        {
        return sqrt(velocity[0] * velocity[0] +
            velocity[1] * velocity[1] +
            velocity[2] * velocity[2]);
        }
        /// Velocity along radial direction, check if returns nan then
        Double_t RadialVelocity() const
        {
        Double_t r = Radius();
        Double_t rv= (velocity[0] * position[0] / r +
            velocity[1] * position[1] / r +
            velocity[2] * position[2] / r);
        if (std::isnan(rv)) return AbsoluteVelocity();
        else return rv;
        }
        /// Velocity tangential to radial direction check if returns nan
        Double_t CircularVelocity() const
        {
            Double_t r = Radius();
            Double_t x = position[1] * velocity[2] - position[2] * velocity[1];
            Double_t y = position[2] * velocity[0] - position[0] * velocity[2];
            Double_t z = position[0] * velocity[1] - position[1] * velocity[0];
            Double_t cv =sqrt(x * x + y * y + z * z) / r;
            if (std::isnan(cv)) return 0.0;
            else return cv;
        }
        /// Radial velocity in cylindrical coordinates.
        Double_t CylRadialVelocity() const
        {
            Double_t p = Phi();
            return velocity[0] * cos(p) + velocity[1] * sin(p);
        }
        /// Tangential velocity in cylindrical coordinates.
        Double_t CylTangentialVelocity() const
        {
            Double_t p = Phi();
            return -velocity[0] * sin(p) + velocity[1] * cos(p);
        }
        /// Angular momentum, J = r * Vt
        Double_t AngularMomentum() const
        {
            return CircularVelocity() * Radius();
        }
        //@}


#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
#ifndef NOMASS
			ar & mass;
#endif
            ar & position;
            ar & velocity;
            ar & id;
            ar & type;
            ar & rho;
        }
        //@}
#endif

    };

    /*
        End of Particle Class
    */


/*!
    \class NBody::GasParticle
    \brief An sph gas particle

    The class is a subclass of \ref NBody::Particle and has several extra quantities commonly need or used in hydro simulations
    These are temperature, log of entropy, pressure, electron or ionization fraction, atomic hydrogen, metallicity, star foramtion rate, internal energy. \n
    This class is in a state of flux as more properties could be added.

*/
    class GasParticle : public Particle
    {
        private:
        /// \name temperature, log of entropy, pressure, electron or ionization fraction, atomic hydrogen, metallicity, star foramtion rate, internal energy what else
        //@{
        Double_t temp, lgS, P, Ne, Nh0, metal, sfr, U;
        //@}
        public:

        /// \name Constructors & Destructors
        //@{
        GasParticle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0,
                Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0,Double_t NH0=0,Double_t Zi=0, Double_t SFR=0,Double_t LGS=0);
        GasParticle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0,
                    Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0,Double_t NH0=0,Double_t Zi=0, Double_t SFR=0,Double_t LGS=0);
        GasParticle(const GasParticle &p);
        //@}

        /// \name Overloaded operators
        //@{
        GasParticle &operator=(const GasParticle &part);
        bool operator==(const GasParticle &p) {return ((Particle::operator==((Particle)p))&&temp==p.temp&&lgS==p.lgS&&P==p.P&&Ne==p.Ne&&Nh0==p.Nh0&&metal==p.metal&&sfr==p.sfr&&U==p.U);}
        bool operator!=(const GasParticle &p) const {return !((Particle::operator==((Particle)p))&&temp==p.temp&&lgS==p.lgS&&P==p.P&&Ne==p.Ne&&Nh0==p.Nh0&&metal==p.metal&&sfr==p.sfr&&U==p.U);}
        //@}

        /// \name Get & Set methods
        //@{
        Double_t GetTemp() const { return temp; }
        void SetTemp(Double_t Temp) { temp = Temp; }
        Double_t GetPressure() const { return P; }
        void SetPressure(Double_t Pi) { P = Pi; }
        Double_t GetU() const { return U; }
        void SetU(Double_t Ui) { U = Ui; }
        Double_t GetNe() const { return Ne; }
        void SetNe(Double_t NE) { Ne = NE; }
        Double_t GetNh0() const { return Nh0; }
        void SetNh0(Double_t NH0) { Nh0 = NH0; }
        Double_t GetZ() const { return metal; }
        void SetZ(Double_t Zi) { metal = Zi; }
        Double_t Getsfr() const { return sfr; }
        void Setsfr(Double_t SFR) { sfr = SFR; }
        Double_t GetEntropy() const { return lgS; }
        void SetEntropy(Double_t lgs) { lgS = lgs; }
        //@}

#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Particle>(*this);
            ar & temp;
            ar & lgS;
            ar & P;
            ar & Ne;
            ar & Nh0;
            ar & metal;
            ar & sfr;
            ar & U;
        }
        //@}
#endif

    };

/*!
    \class NBody::StarParticle
    \brief An sph star particle

    The class is a subclass of \ref NBody::Particle and has several extra quantities commonly need or used in hydro simulations
    These are formation time and metallicity. \n
    This class is in a state of flux as more properties could be added.

*/
    class StarParticle : public Particle
    {
        private:

        ///formation time
        Double_t tform;
        ///metallicity
        Double_t metal;

        public:
        /// \name Constructors & Destructors
        //@{
        StarParticle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0);
        StarParticle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, PARTIDTYPE ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0);
        StarParticle(const StarParticle &p);
        //@}

        /// \name Overloaded operators
        //@{
        StarParticle &operator=(const StarParticle &part);
        bool operator==(const StarParticle &p) {return ((Particle::operator==((Particle)p))&&tform==p.tform&&metal==p.metal);}
        bool operator!=(const StarParticle &p) const { return !((Particle::operator==((Particle)p))&&tform==p.tform&&metal==p.metal);}
        //@}

        /// \name Get & Set methods
        //@{
        Double_t GetFormationTime() const { return tform; }
        void SetFormationTime(Double_t Tform) { tform = Tform; }
        Double_t GetZ() const { return metal; }
        void SetZ(Double_t zz) { metal = zz; }
        //@}

#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Particle>(*this);
            ar & tform;
            ar & Z;
        }
        //@}
#endif

    };
    /*
        End of Particle SubClasses
    */

    ///\name useful comparison functions for qsort involving particles
    //@{
    ///sort in ascending particle pid
    int PIDCompare (const void *a, const void *b);
    ///sort in ascending particle id
    int IDCompare (const void *a, const void *b);
    ///sort in ascending particle radius
    int RadCompare (const void *a, const void *b);
    ///sort in ascending particle type
    int TypeCompare (const void *a, const void *b);
    ///sort in ascending particle potential
    int PotCompare (const void *a, const void *b);
    ///sort in ascending particle density
    int DenCompare (const void *a, const void *b);

    ///sort in ascending particle pid for std::sort vector inferface
    bool PIDCompareVec (const Particle &a, const Particle &b);
    ///sort in ascending particle id
    bool IDCompareVec (const Particle &a, const Particle &b);
    ///sort in ascending particle radius
    bool RadCompareVec (const Particle &a, const Particle &b);
    ///sort in ascending particle type
    bool TypeCompareVec (const Particle &a, const Particle &b);
    ///sort in ascending particle potential
    bool PotCompareVec (const Particle &a, const Particle &b);
    ///sort in ascending particle density
    bool DenCompareVec (const Particle &a, const Particle &b);
    //@}
}

#ifdef USEBOOSTMPI
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::Particle> : mpl::true_ { };
        } }
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::GasParticle> : mpl::true_ { };
        } }
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::StarParticle> : mpl::true_ { };
        } }
#endif

#endif // PARTICLE_H
