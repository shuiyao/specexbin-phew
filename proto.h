#ifdef PHEW
void set_random_numbers(void);
double get_random_number(unsigned int id);
#endif

int RotateCoords(double, double, double, double, double);
int InverseRotateCoords(double, double, double, double, double);
int binarysearch(double, double *, int);

/* initions */
//int InitIons()

/* file_io */
int Check_Z_File();
int Open_Snapshot(char *);

/* read_hdf5 */
int readhdf5(char*, int);
double readunits(char*, int, char*);

/* contsmoothspec */
//int ContSmoothSpec();

/* tau */
//int Tau();

/* outtau */
//int OutTau();
