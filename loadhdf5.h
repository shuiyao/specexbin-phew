#ifndef _LOADHDF5_H
#define _LOADHDF5_H

#define MAX_LEN_FILENAME 500
#define MAXDIM 3
#define NMETALSHDF5 11

#define UNIT_L 3.085678e21
#define UNIT_V 1.e5
#define UNIT_M 1.989e43

/* Variables */
extern double unit_Time;
extern double unit_Length;
extern double unit_Density;
extern double unit_Mass;
extern double unit_Velocity;
extern double unit_Temp;

struct gadget_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 

  int      flag_stellarage;
  int      flag_metals;

  unsigned int      npartTotalHighWord[6];

  int      flag_entropy_instead_u;
  int      flag_doubleprecision;

  int flag_potential;
  int flag_fH2;

  int flag_tmax;
  int flag_delaytime;

  int flag_lpt_ics;
  float flag_lpt_scalingcator;

  char     fill[32];  /* fills to 256 Bytes */
} ;

struct particle_data
{
  float  Pos[3];
  float  Vel[3];
  int ID;
  float Mass, Rho, Ne, Nh, Hsml, metal[NMETALSHDF5], fH2, Tmax;
  float Sfr; // StellarFormationTime for star particles
  float Temp;
  float DelayTime;
  int    Flag;
#ifdef PHEW
  int Key;
  float WindMass;
  float Mcloud;
  float Rcloud;
  float LastSFTime;
  float Vinit;
#endif  
};

extern struct particle_data *P;
extern struct gadget_header gheader;

/* Functions */
int load_hdf5(char *snap, int itype);
int allocate_memory(void);
void tipsyunits();

#endif // _LOADHDF5_H
