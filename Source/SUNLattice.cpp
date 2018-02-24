#include "SUNLattice.h"

namespace SUN
{	
 
matrix3* Lattice;
int acc(0);
int dec(0);
vector<LatticePosition> Basis;

SUNLattice::SUNLattice(LatticeSettings& Settings):settings(Settings)
{
	Basis.push_back(LatticePosition(1,0,0,0,settings));
	Basis.push_back(LatticePosition(0,1,0,0,settings));
	Basis.push_back(LatticePosition(0,0,1,0,settings));
	Basis.push_back(LatticePosition(0,0,0,1,settings));
	random= new Random(settings);
	Lattice = new matrix3[settings.Volume()*4];
	if(settings.Condition == StartCondition::Cold)
	{
		//#pragma omp parallel for
		for (int i = 0; i < settings.Volume()*4 ; i++)
		{
			Lattice[i]= matrix3(MatrixInitType::Unit);
		}
		
	}
	else
	{
		//#pragma omp parallel for
		for (int i = 0; i < settings.Volume()*4 ; i++)
		{
			Lattice[i]= random->GetUnitaryMatrix3();
		}
		
	}
}


int SUNLattice::CoordinatesToIndex(int t, int x, int y, int z, int mu)
{
	return 	mu+4*t\
			+4*settings.t_LattSize*x\
			+4*settings.t_LattSize*settings.s_LattSize*y\
			+4*settings.t_LattSize*settings.s_LattSize*settings.s_LattSize*z;
}
int SUNLattice::LatticePositionToIndex(LatticePosition pos,int mu)
{
	return SUNLattice::CoordinatesToIndex(pos.Pos[0],pos.Pos[1],pos.Pos[2],pos.Pos[3],mu);
}
LatticePosition SUNLattice::IndexToLatticePosition(int n, int * mu)
{
	int p[4];
	p[3]=n/(4*settings.t_LattSize*settings.s_LattSize*settings.s_LattSize);
	n=n%(4*settings.t_LattSize*settings.s_LattSize*settings.s_LattSize);
	p[2]=n/(4*settings.t_LattSize*settings.s_LattSize);
	n=n%(4*settings.t_LattSize*settings.s_LattSize);
	p[1]=n/(4*settings.t_LattSize);
	n=n%(4*settings.t_LattSize);
	p[0]=n/(4);
	n=n%4;
	*mu = n;
	return 	LatticePosition(p[0],p[1],p[2],p[3],settings);
}

SUNLattice::~SUNLattice()
{
	delete[] Lattice;
}

matrix3 SUNLattice::GetLink(LatticePosition pos, int mu)
{
	return Lattice[SUNLattice::LatticePositionToIndex(pos,mu)];
}

matrix3	SUNLattice::Staple(LatticePosition n, int mu)
{
	matrix3 A= matrix3(MatrixInitType::Zero);
	for (int nu = 0; nu < 4; nu++)
	{
		if(nu==mu)
			continue;
		LatticePosition pos = n+Basis[mu];
		A=A	+GetLink(n+Basis[mu],nu)*(GetLink(n+Basis[nu],mu).dagger())*(GetLink(n,nu).dagger()) \
			+(GetLink(pos-Basis[nu],nu).dagger())*(GetLink(n-Basis[nu],mu).dagger())*GetLink(n-Basis[nu],nu);
		
	}
	return A;
	
}

double SUNLattice::LocalSChange(LatticePosition n, int mu, MATRIX& trial)
{
	return real(((GetLink(n,mu)-trial)*Staple(n,mu)).trace())*settings.beta/GROUPDIM;
}

void SUNLattice::MonteCarloStep()
{
	LatticePosition pos = random->GetLatticePosition();
	int mu = random->GetInt(4);
	MATRIX trial = random->GetUnitaryMatrix3();
	if(random->GetDouble()<exp(-LocalSChange(pos,mu,trial)))
	{
		acc++;
		Lattice[LatticePositionToIndex(pos,mu)]=trial;
	}
	else
		dec++;
}

void SUNLattice::MonteCarloUpdate()
{
	acc=0;
	dec=0;
	//#pragma omp parallel for
	//#pragma omp parallel for
	for (int i = 0; i < settings.Volume()*4; i++)
	{
		SUNLattice::MonteCarloStep();
	}
}

MATRIX SUNLattice::Plaquette(LatticePosition n, int mu, int nu, int a, int b)
{
	matrix3 xSide(MatrixInitType::Unit);
	matrix3 ySide(MatrixInitType::Unit);
	matrix3 xSideBack(MatrixInitType::Unit);
	matrix3 ySideBack(MatrixInitType::Unit);
	for (int i = 0; i < a; i++)
	{
		xSide = xSide * GetLink(n+Basis[mu]*i,mu);
		xSideBack = xSideBack*(GetLink(n+Basis[mu]*(a-i-1)+Basis[nu]*b,mu).dagger());
	}
	
	for (int i = 0; i < b; i++)
	{
		ySide = ySide*GetLink(n+Basis[mu]*a+Basis[nu]*i,nu);
		ySideBack = ySideBack*(GetLink(n+Basis[nu]*(b-i-1),nu).dagger());
	}
	
	
	return xSide*ySide*xSideBack*ySideBack;
}

double SUNLattice::WilsonLoop(int a, int b)
{
	double res = 0.0;
	int count = 0;
	//#pragma omp parallel for reduction(+ : res , count) collapse(4)
	forlatt
	{
		LatticePosition pos(t,x,y,z,settings);
		for (int mu = 0; mu < 4; mu++)
		{
			for (int nu = 0; nu < 4; nu++)
			{
				if(nu==mu)
					continue;
				count++;
				res+=real(Plaquette(pos,mu,nu,a,b).trace());
			}					
		}		
	}
	return res /(GROUPDIM*count);
}

Complex SUNLattice::PlaqPlaqCorrelation(int tdiff)
{
	int size = 1;
	double im, re;
	im = 0.0;
	re = 0.0;
	Complex res = 0.0;
	int count = 0;
	#pragma omp parallel for reduction(+ : im, re , count) collapse(4)
	forlatt
	{
		LatticePosition pos(t,x,y,z,settings);
		for (int mu = 1; mu < 4; mu++)
		{
			for (int nu = 1; nu < 4; nu++)
			{
				if(nu==mu)
					continue;
				count++;
				res=Plaquette(pos,mu,nu,size,size).trace()*Plaquette(pos+Basis[0]*tdiff,mu,nu,size,size).trace();
				im += imag(res);
				re += real(res);
			}					
		}		
	}
	return Complex(re,im) /(GROUPDIM*GROUPDIM*count);
}

Complex SUNLattice::PolyakovLoop(int x, int y, int z)
{
	matrix3 Loop(MatrixInitType::Unit);
	for (int t = 0; t < settings.t_LattSize; t++)
	{
		LatticePosition pos(t,x,y,z,settings);
		Loop = Loop * GetLink(pos, 0);
	}
	return Loop.trace();
	
}

int SUNLattice::PolyakovPosToIndex(int x, int y, int z)
{
	return (x%settings.s_LattSize)+(y%settings.s_LattSize)*settings.s_LattSize+(z%settings.s_LattSize)*settings.s_LattSize*settings.s_LattSize;
}

void SUNLattice::PolyakovIndexToPos(int index, int &x, int &y, int &z)
{
	z = index /(settings.s_LattSize*settings.s_LattSize);
	index = index %(settings.s_LattSize*settings.s_LattSize);
	y = index /settings.s_LattSize;
	x = index % settings.s_LattSize;
}

void SUNLattice::PolyakovLoopCorrelation()
{
	int avgcound = 0;
//	int maxSeperation = settings.s_LattSize/2;
	double Result ;
	Complex* PolyakovLoops = new Complex[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize];
	for (int z = 0; z < settings.s_LattSize; z++)
	{
		for (int y = 0; y < settings.s_LattSize; y++)
		{
			for (int x = 0; x < settings.s_LattSize; x++)
			{
				PolyakovLoops[SUNLattice::PolyakovPosToIndex(x,y,z)]= SUNLattice::PolyakovLoop(x,y,z);
			}
		}
	}

		Result=0.0;
		for (int z = 0; z < settings.s_LattSize; z++)
		{
			for (int y = 0; y < settings.s_LattSize; y++)
			{
				for (int x = 0; x < settings.s_LattSize; x++)
				{
					Result+=real(PolyakovLoops[SUNLattice::PolyakovPosToIndex(x,y,z)]*conj(PolyakovLoops[SUNLattice::PolyakovPosToIndex(x+2,y,z)]));	
					avgcound ++;				
				}
			}
		}
		Result/=avgcound;
		printf("%f	",Result);
	printf("\n");
	delete[] PolyakovLoops;
	return;
}

void SUNLattice::PolyakovRun()
{
	int thermN = 1000;
	int N = 100;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		printf("Thermalising : %i / %i\n", i+1, thermN);
	}
	for (int i = 0; i < N; i++)
	{
		SUNLattice::MonteCarloUpdate();
		SUNLattice::PolyakovLoopCorrelation();
	}
	
}

void SUNLattice::AvgActionPerPlaqueRun()
{
}

void SUNLattice::WilsonLoopRun(const char * filename)
{
	double res;
	FILE * f = fopen(filename,"a");
	fprintf(f,"## Beta : %f	L : %i	T : %i\n", settings.beta,settings.s_LattSize,settings.t_LattSize);
	fclose(f);
	
	for (int i = 0; i < 1000; i++)
	{
		SUNLattice::MonteCarloUpdate();
		
		//printf("Thermalising : %i / %i\n", i+1, 1000);
		
	}
	
	for (int i = 0; i < 100000; i++)
	{
		//if(i%1000==0)
		//	printf("%i\n",i);
		SUNLattice::MonteCarloUpdate();
		f = fopen(filename,"a");
		res=WilsonLoop(1,1);
		fprintf(f,"%i	%f\n",i,res);
		fclose(f);
	}
	return;
	
	
	int thermN = 0;
	int N = 2000;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		printf("Thermalising : %i / %i\n", i+1, thermN);
	}
	for (int i = 0; i < N; i++)
	{
		SUNLattice::MonteCarloUpdate();
		printf("%f	%f\n",SUNLattice::WilsonLoop(1,1),SUNLattice::WilsonLoop(2,2));
	}
}

void SUNLattice::GenQuarkPot()
{
	int maxsize = settings.s_LattSize;
	int thermN = 2000;
	int N = 2000;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		fprintf( stderr,"#Thermalising : %i / %i\n", i+1, thermN);
	}
	fprintf( stdout,"#");
	for (int x = 1; x < maxsize; x++)
			fprintf( stdout,"%i	",x);
	fprintf( stdout,"\n");
	
	for (int i = 0; i < N; i++)
	{
		//fprintf( stdout,"%i	",i);
		SUNLattice::MonteCarloUpdate();
		for (int x = 1; x < maxsize; x++)
				fprintf( stdout,"%f	",SUNLattice::WilsonLoop(1,x));
		
		
		fprintf( stdout,"\n");
	}
}

void SUNLattice::PlaqPlaqRun()
{
	int thermN = 2000;
	int N = 20000000;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		fprintf( stderr,"#Thermalising : %i / %i\n", i+1, thermN);
	}
	Complex res;
	for (int i = 0; i < N; i++)
	{
		printf("%i",i);
		SUNLattice::MonteCarloUpdate();
		for (int sep = 1; sep < 15; sep++)
		{
			res = PlaqPlaqCorrelation(sep);
			printf("	%f	%f",real(res),imag(res));
		}
		
		printf("\n");
	}
}
/*
void SUNLattice::SaveLattice()
{
}
void SUNLattice::LoadLattice()
{
}

const int rAvg 	= 2;
const int n_up0 = 10;
const int n_ms0 = 10;
//const int n_up1 = 10;
//const int n_ms1 = 10;
TwoLinkOperator * opAvgs0;
TwoLinkOperator * opAvgs1;
void SUNLattice::TimeSiliceUpdate(int x0, int y0)
{
	int tdiff = y0-x0;
	
	int mu;
	for (int i = 0; i < settings.Volume()/settings.t_LattSize*tdiff; i++)
	{
		LatticePosition pos = random->GetLatticePosition();
		pos.Pos[0]=x0+random->GetInt(tdiff+1);
		mu = random->GetInt(4);
		if(pos.Pos[0]==x0 || pos.Pos[0]==y0) 	//Lattice gauge theory on a time-slice can be studied independently of the surrounding
		{
			if(mu!=0)
			{
				i--;
				continue;
				}//lattice if the spatial link variables at the boundaries are held fixed
		}
		//printf("%i	%i	%i	%i	%i\n",pos.Pos[0],pos.Pos[1],pos.Pos[2],pos.Pos[3],mu);
		MATRIX trial = random->GetUnitaryMatrix3();
		if(random->GetDouble()<exp(-LocalSChange(pos,mu,trial)))
		{
			acc++;
			Lattice[LatticePositionToIndex(pos,mu)]=trial;
		}
		else
			dec++;
		}
}

int SUNLattice::ThreeCoordinatesToIndex(int x, int y, int z, int timeSlice, int direction, int SliceCount)
{
	return (x%settings.s_LattSize)\
			+(y%settings.s_LattSize)*settings.s_LattSize\
			+(z%settings.s_LattSize)*settings.s_LattSize*settings.s_LattSize\
			+(timeSlice%SliceCount)*settings.s_LattSize*settings.s_LattSize*settings.s_LattSize\
			+(direction-1)*SliceCount*settings.s_LattSize*settings.s_LattSize*settings.s_LattSize ;
}

void SUNLattice::AvgLevel(int level)
{
	int r = 3;
	
	int TimeSliceCount=0;
	if(level ==0)
	{
		TimeSliceCount = settings.s_LattSize/2;
		opAvgs0=new TwoLinkOperator[TimeSliceCount*settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		for (int i = 0; i < TimeSliceCount; i++)
		{
			//printf("Avging  ts %i / %i in level 0\n",i+1,TimeSliceCount);
			for (int j = 0; j < n_ms0; j++)
			{
				for (int x = 0; x < settings.s_LattSize; x++)
				{
					for (int y = 0; y < settings.s_LattSize; y++)
					{
						for (int z = 0; z < settings.s_LattSize; z++)
						{
							for (int  dir = 1; dir < 4; dir++)
							{
								LatticePosition pos(2*i,x,y,z);
								opAvgs0[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)]=opAvgs0[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)]+TwoLinkOperator(GetLink(pos,0),GetLink(pos+Basis[dir]*r,0))*TwoLinkOperator(GetLink(pos+Basis[0],0),GetLink(pos+Basis[0]+Basis[dir]*r,0));
							}
						}
					}
				}
				for (int up = 0; up< n_up0; up++)
					SUNLattice::TimeSiliceUpdate(i*2, i*2+2);
			}
			for (int x = 0; x < settings.s_LattSize; x++)
			{
				for (int y = 0; y < settings.s_LattSize; y++)
				{
					for (int z = 0; z < settings.s_LattSize; z++)
					{
						for (int  dir = 1; dir < 4; dir++)
						{
							opAvgs0[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)].DivideByScalar((double)n_ms0);
						}
					}
				}
			}
		}
	}
	if(level ==1)
	{
		TimeSliceCount = settings.s_LattSize/4;
		
		opAvgs1=new TwoLinkOperator[TimeSliceCount*settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		for (int j = 0; j < n_ms1; j++)
		{
			printf("#Avging %i / %i in level 1\n\n",j+1,n_ms1);
			for (int up = 0; up< n_up1; up++)
				AvgLevel(0);
			for (int i = 0; i < TimeSliceCount; i++)
			{
				for (int x = 0; x < settings.s_LattSize; x++)
				{
					for (int y = 0; y < settings.s_LattSize; y++)
					{
						for (int z = 0; z < settings.s_LattSize; z++)
						{
							for (int  dir = 1; dir < 4; dir++)
							{
								opAvgs1[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)]=opAvgs1[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)]+opAvgs0[ThreeCoordinatesToIndex(x,y,z,i*2,dir,TimeSliceCount*2)]*opAvgs0[ThreeCoordinatesToIndex(x,y,z,i*2+1,dir,TimeSliceCount*2)];
							}
						}
					}
				}
			}
			delete[] opAvgs0;
		}
		for (int i = 0; i < TimeSliceCount; i++)
		{
			for (int x = 0; x < settings.s_LattSize; x++)
			{
				for (int y = 0; y < settings.s_LattSize; y++)
				{
					for (int z = 0; z < settings.s_LattSize; z++)
					{
						for (int  dir = 1; dir < 4; dir++)
						{
							opAvgs1[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)].DivideByScalar((double)n_ms1);
						}
					}
				}
			}
		}
	}
		
}

void SUNLattice::MultiRun()
{
	int TimeSliceCount = settings.s_LattSize/4;
	int thermN = 1000;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		fprintf( stderr,"#Thermalising : %i / %i\n", i+1, thermN);
	}
	//int r = 2;
	double traceRes;
	//TwoLinkOperator * opAvgs;
	TwoLinkOperator res(matrix3(MatrixInitType::Unit),matrix3(MatrixInitType::Unit));
	for(int m =0; m< 100;m++)
	{
		SUNLattice::AvgLevel(1);
		traceRes = 0.0;
		for (int x = 0; x < settings.s_LattSize; x++)
			{
				for (int y = 0; y < settings.s_LattSize; y++)
				{
					for (int z = 0; z < settings.s_LattSize; z++)
					{
						for (int  dir = 1; dir < 4; dir++)
						{
							TwoLinkOperator res(matrix3(MatrixInitType::Unit),matrix3(MatrixInitType::Unit));
							for (int i = 0; i < TimeSliceCount; i++)
							{
								res = res * opAvgs1[ThreeCoordinatesToIndex(x,y,z,i,dir,TimeSliceCount)];
							}
							traceRes+= real(res.trace());
						}	
					}
				}
			}
		SUNLattice::MonteCarloUpdate();
		SUNLattice::MonteCarloUpdate();
		printf("%e	\n",traceRes/pow(settings.s_LattSize,3)/3);
		delete[] opAvgs1;
	}
	
}

const int n_up1 = 5;
const int n_ms1 = 5;
const int n_up2 = 5;
const int n_ms2 = 5;
const int n_up4 = 5;
const int n_ms4 = 5;
const int RAvg 	= 2;


int SUNLattice::ThreeCoordinatesToIndex(int x, int y, int z, int direction)
{
	return 	(x%settings.s_LattSize)\
			+(y%settings.s_LattSize)*settings.s_LattSize\
			+(z%settings.s_LattSize)*settings.s_LattSize*settings.s_LattSize\
			+(direction-1)*settings.s_LattSize*settings.s_LattSize*settings.s_LattSize ;
}


void SUNLattice::Avg1aSlices(int slice, TwoLinkOperator * result)
{
	for (int measure = 0; measure < n_ms1; measure++)
	{	
		for (int i = 0; i < n_up1; i++)
		{
			SUNLattice::TimeSiliceUpdate(slice,slice+1);
		}
		for (int x = 0; x < settings.s_LattSize; x++)
		{
			for (int y = 0; y < settings.s_LattSize; y++)
			{
				for (int z = 0; z < settings.s_LattSize; z++)
				{
					for (int direction = 1; direction <4 ; direction++)
					{
						LatticePosition pos(slice,x,y,z,settings);
						result[ThreeCoordinatesToIndex(x,y,z,direction)]=result[ThreeCoordinatesToIndex(x,y,z,direction)]+TwoLinkOperator(GetLink(pos,0),GetLink(pos+Basis[direction]*RAvg,0));
					}
				}
			}
		}
	}
	for (int x = 0; x < settings.s_LattSize; x++)
	{
		for (int y = 0; y < settings.s_LattSize; y++)
		{
			for (int z = 0; z < settings.s_LattSize; z++)
			{
				for (int direction = 1; direction <4 ; direction++)
				{
					result[ThreeCoordinatesToIndex(x,y,z,direction)].DivideByScalar((double) n_ms1);
				}
			}
		}
	}
}

void SUNLattice::Avg2aSlices(int slice, TwoLinkOperator * result)
{
	//printf("# Slice %i / %i\n",slice+1, settings.s_LattSize/2);
	TwoLinkOperator * slice1;
	TwoLinkOperator * slice2;
	
	for (int measure = 0; measure < n_ms2; measure++)
	{	
		for (int i = 0; i < n_up2; i++)
		{
			SUNLattice::TimeSiliceUpdate(slice,slice+2);
		}
		slice1 = new TwoLinkOperator[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		slice2 = new TwoLinkOperator[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		SUNLattice::Avg1aSlices(slice, slice1);
		SUNLattice::Avg1aSlices(slice+1, slice2);
		for (int x = 0; x < settings.s_LattSize; x++)
		{
			for (int y = 0; y < settings.s_LattSize; y++)
			{
				for (int z = 0; z < settings.s_LattSize; z++)
				{
					for (int direction = 1; direction <4 ; direction++)
					{
						result[ThreeCoordinatesToIndex(x,y,z,direction)]=result[ThreeCoordinatesToIndex(x,y,z,direction)]+slice1[ThreeCoordinatesToIndex(x,y,z,direction)]*slice2[ThreeCoordinatesToIndex(x,y,z,direction)];
					}
				}
			}
		}
		delete[] slice1;
		delete[] slice2;
		
	}
	for (int x = 0; x < settings.s_LattSize; x++)
	{
		for (int y = 0; y < settings.s_LattSize; y++)
		{
			for (int z = 0; z < settings.s_LattSize; z++)
			{
				for (int direction = 1; direction <4 ; direction++)
				{
					result[ThreeCoordinatesToIndex(x,y,z,direction)].DivideByScalar((double) n_ms2);
				}
			}
		}
	}		
}

void SUNLattice::Avg4aSlices(int slice, TwoLinkOperator * result)
{
	TwoLinkOperator * slice1;
	TwoLinkOperator * slice2;
	
	for (int measure = 0; measure < n_ms4; measure++)
	{	
		for (int i = 0; i < n_up4; i++)
		{
			SUNLattice::TimeSiliceUpdate(slice,slice+4);
		}
		slice1 = new TwoLinkOperator[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		slice2 = new TwoLinkOperator[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		SUNLattice::Avg2aSlices(slice, slice1);
		SUNLattice::Avg2aSlices(slice+2, slice2);
		for (int x = 0; x < settings.s_LattSize; x++)
		{
			for (int y = 0; y < settings.s_LattSize; y++)
			{
				for (int z = 0; z < settings.s_LattSize; z++)
				{
					for (int direction = 1; direction <4 ; direction++)
					{
						result[ThreeCoordinatesToIndex(x,y,z,direction)]=result[ThreeCoordinatesToIndex(x,y,z,direction)]+slice1[ThreeCoordinatesToIndex(x,y,z,direction)]*slice2[ThreeCoordinatesToIndex(x,y,z,direction)];
					}
				}
			}
		}
		delete[] slice1;
		delete[] slice2;
		
	}
	for (int x = 0; x < settings.s_LattSize; x++)
	{
		for (int y = 0; y < settings.s_LattSize; y++)
		{
			for (int z = 0; z < settings.s_LattSize; z++)
			{
				for (int direction = 1; direction <4 ; direction++)
				{
					result[ThreeCoordinatesToIndex(x,y,z,direction)].DivideByScalar((double) n_ms4);
				}
			}
		}
	}		
}

double SUNLattice::TraceAvgComlete()
{
	TwoLinkOperator * sliceAverages;
	vector<TwoLinkOperator> result(settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3,TwoLinkOperator(matrix3(MatrixInitType::Unit),matrix3(MatrixInitType::Unit)));
	for (int slice = 0; slice < settings.s_LattSize/4; slice++)
	{
		//printf("# Slice %i / %i\n",slice+1, settings.s_LattSize/4);
		sliceAverages = new TwoLinkOperator[settings.s_LattSize*settings.s_LattSize*settings.s_LattSize*3];
		Avg4aSlices(4*slice,sliceAverages);
		for (int x = 0; x < settings.s_LattSize; x++)
		{
			for (int y = 0; y < settings.s_LattSize; y++)
			{
				for (int z = 0; z < settings.s_LattSize; z++)
				{
					for (int direction = 1; direction <4 ; direction++)
					{
						result[ThreeCoordinatesToIndex(x,y,z,direction)]= result[ThreeCoordinatesToIndex(x,y,z,direction)]*sliceAverages[ThreeCoordinatesToIndex(x,y,z,direction)];
					}
				}
			}
		}
		delete[] sliceAverages;
	} 
	double res = 0.0;
	for (int x = 0; x < settings.s_LattSize; x++)
	{
		for (int y = 0; y < settings.s_LattSize; y++)
		{
			for (int z = 0; z < settings.s_LattSize; z++)
			{
				for (int direction = 1; direction <4 ; direction++)
				{
					res+=real(result[ThreeCoordinatesToIndex(x,y,z,direction)].trace());
				}
			}
		}
	}
	return res / pow(settings.s_LattSize,3)/3;
}
void SUNLattice::MultiRun2()
{
	int thermN = 1000;
	for (int i = 0; i < thermN; i++)
	{
		SUNLattice::MonteCarloUpdate();
		fprintf( stderr,"#Thermalising : %i / %i\n", i+1, thermN);
	}
	for (int i = 0; i < 2000; i++)
	{
		for (int mc = 0; mc < 100; mc++)
		{
			SUNLattice::MonteCarloUpdate();
		}
		printf("%e\n",SUNLattice::TraceAvgComlete());
	}
	
	
}
*/
void SUNLattice::Run(const char * filename)
{
	//TimeSiliceUpdate(0,2);
	WilsonLoopRun(filename);
	//printf("%f",real(SUNLattice::PolyakovLoop(1,1,1)));
	return;
}
}
