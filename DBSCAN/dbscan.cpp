#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unordered_map>
#include <vector>
#include <stack>
#include <algorithm>
using std::unordered_map;
using std::vector;
using std::stack;

// const parameter
double Pi = 3.1415926;
double eps = 1e-10;
int scale = 1e6;

// self define structure
struct Particle
{
	double *X, *Y, *Z;
	int *ID, *Flag, *inClst; // Flag means visit flag; inClst is the cluster id of the particle
	int Nr;
	int *Core, *Border, *Noise;
};

struct Cluster
{
	vector<int> idList;
	double diam;
	int size;
	double CooNum;
};

// read from para
double Skin;
int MinPts;
int outputVTK;
double radius;
double L[3][2];
double FileTime;

// global parameter
int startStep, endStep, nowStep, dt, filenum, debugStep, continueStep;
double dp, dp_skin_2, Rskin;
double PrsL[3], HalfPrsL[3], Grid[3], T[3], r_Grid[3];
int N[3], B[3];
int vect[27][3], VNum = 27;

int *bin, *binlen, BLen, BNum;
int *binxid, *binyid, *binzid;
int *NBList, *NBLength, NLen;
int *CooNum;
double *NBDis;
double *ClstDis;
int ClusterNum, *ClusterLen;
vector<Cluster*> dbClst;
Particle p;

FILE* c_debug;

void ReadDBSCANInfo()
{
	FILE *fp;
	char sN[200] = "dbscanpara.ini";

	fp = fopen(sN, "r");
	fscanf(fp, "radius=%lf\n", &radius);
	fscanf(fp, "outputVTK=%d\n", &outputVTK);
	fscanf(fp, "Skin=%lf\n", &Skin);
	fscanf(fp, "MinPts=%d\n", &MinPts);
	fscanf(fp, "box0=%lf\t%lf\t%lf\n", &L[0][0], &L[1][0], &L[2][0]);
	fscanf(fp, "box1=%lf\t%lf\t%lf\n", &L[0][1], &L[1][1], &L[2][1]);
	fscanf(fp, "FileTime=%lf\n", &FileTime);

	fclose(fp);

	dp = 2.0*radius;
	printf("Read para done\n");
}

void InitializeNeighbourList()
{
	Rskin = Skin*dp;
	dp_skin_2 = (dp + Rskin)*(dp + Rskin);

	double G = 5 * dp + Rskin;
	for (int i = 0; i<3; i++)
	{
		PrsL[i] = (L[i][1] - L[i][0]);
		HalfPrsL[i] = PrsL[i] * 0.5;

		int binnum = int(PrsL[i] / G + eps);
		if (binnum == 0)
			binnum = 1;

		Grid[i] = PrsL[i] / double(binnum);
		r_Grid[i] = 1.0 / Grid[i];

		T[i] = PrsL[i];
		B[i] = binnum;
	}

	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			for (int k = 0; k<3; k++)
			{
				int tmp = i * 9 + j * 3 + k;
				vect[tmp][0] = j - 1;
				vect[tmp][1] = k - 1;
				vect[tmp][2] = i - 1;
			}
	printf("Init NBL done\n");
}

void InitializeParticleArray(Particle &a, int Len)
{
	a.X = (double *)malloc(sizeof(double)*Len);
	a.Y = (double *)malloc(sizeof(double)*Len);
	a.Z = (double *)malloc(sizeof(double)*Len);
	a.ID = (int *)malloc(sizeof(int)*Len);
	a.Flag = (int *)malloc(sizeof(int)*Len);
	a.inClst = (int *)malloc(sizeof(int)*Len);
	a.Core = (int *)malloc(sizeof(int)*Len);
	a.Border = (int *)malloc(sizeof(int)*Len);
	a.Noise = (int *)malloc(sizeof(int)*Len);

	if (!a.X || !a.Y || !a.Z || !a.Flag || !a.ID || !a.inClst || !a.Border || !a.Noise || !a.Core)
	{
		printf("InitializeParticleArray error, failed to allocate memory.\n");
		exit(1);
	}

	memset(a.X, 0, Len * sizeof(double));
	memset(a.Y, 0, Len * sizeof(double));
	memset(a.Z, 0, Len * sizeof(double));
	memset(a.ID, 0, Len * sizeof(int));
	memset(a.Flag, 0, Len * sizeof(int));
	memset(a.inClst, 0, Len * sizeof(int));
	memset(a.Border, 0, Len * sizeof(int));
	memset(a.Core, 0, Len * sizeof(int));
	memset(a.Noise, 0, Len * sizeof(int));

	printf("Init patic array done\n");
}

void AllocMem(int NP)
{
	// nblist
	NLen = int(8 * (dp + Rskin)*(dp + Rskin)*(dp + Rskin) / dp / dp / dp + 50 + eps);
	NBLength = (int *)malloc(sizeof(int)*NP);
	NBList = (int *)malloc(sizeof(int)*(NP*NLen));
	NBDis = (double *)malloc(sizeof(double)*(NP*NLen));
	if (!NBLength || !NBList || !NBDis)
	{
		printf("AllocMem error, failed to allocate nblist memory.\n");
		exit(1);
	}
	memset(NBLength, 0, NP * sizeof(int));
	memset(NBList, 0, (NP*NLen) * sizeof(int));
	memset(NBDis, 0, (NP*NLen) * sizeof(double));

	// bin
	binxid = (int *)malloc(sizeof(int)*NP);
	binyid = (int *)malloc(sizeof(int)*NP);
	binzid = (int *)malloc(sizeof(int)*NP);
	if (!binyid || !binxid || !binzid)
	{
		printf("AllocMem error, failed to allocate bin memory.\n");
		exit(1);
	}
	memset(binxid, 0, NP * sizeof(int));
	memset(binxid, 0, NP * sizeof(int));
	memset(binxid, 0, NP * sizeof(int));

	BLen = int(Grid[0] * Grid[1] * Grid[2] * 6 / Pi / dp / dp / dp + 30 + eps);
	BNum = B[0] * B[1] * B[2];
	bin = (int *)malloc(sizeof(int)*(BNum*BLen));
	binlen = (int *)malloc(sizeof(int)*(BNum));
	if (!bin || !binlen)
	{
		printf("InitializeNeighbourList error, failed to allocate bin memory.\n");
		exit(1);
	}
	memset(binlen, 0, BNum * sizeof(int));
	memset(bin, 0, BNum*BLen * sizeof(int));

	// clst
	ClusterLen = (int *)malloc(sizeof(int)*NP);
	CooNum = (int *)malloc(sizeof(int)*NP);
	ClstDis = (double *)malloc(sizeof(double)*NP);
	if (!ClusterLen || !CooNum || !ClstDis)
	{
		printf("AllocMem error, failed to allocate cluster memory.\n");
		exit(1);
	}
	memset(ClusterLen, 0, NP * sizeof(int));
	memset(CooNum, 0, NP * sizeof(int));
	memset(ClstDis, 0, NP * sizeof(double));

	printf("AllocMem done.\n");
}

void ReadVTK(char* name, Particle &patic)
{
	char dump[1000];

	FILE* fp = fopen(name, "r");
	fscanf(fp, "%[^\n]\n", &dump);
	fscanf(fp, "\n");
	fscanf(fp, "%[^\n]\n", &dump);
	fscanf(fp, "%[^\n]\n", &dump);
	fscanf(fp, "%s  %d %s \n", dump, &patic.Nr, dump);

	printf("NP=%d\n", patic.Nr);

	InitializeParticleArray(patic, patic.Nr);

	for (int i = 0; i < patic.Nr; i++)
	{
		fscanf(fp, "%lf\t%lf\t%lf\n", &patic.X[i], &patic.Y[i], &patic.Z[i]);
		//if(i==0) printf("px,py,pz:(%lf,%lf,%lf)\n",px[i],py[i],pz[i]);
		//if(i==NP-1) printf("px,py,pz:(%lf,%lf,%lf)\n",px[i],py[i],pz[i]);
	}

	fclose(fp);
	AllocMem(patic.Nr);
	printf("Read %s done\n", name);
}

int Locate2(int a, int b, int num)
{
	return (b*num + a);
}

int Locate3(int a, int b, int c)
{
	return (c*(B[0] * B[1]) + b*(B[0]) + a);
}

int Locate4(int a, int b, int c, int d)
{
	return (d*(BNum)+c*(B[0] * B[1]) + b*(B[0]) + a);
}

void Updatebin(int start, int Len, Particle &patic)
{
	for (int i = start; i < start + Len; i++)
	{
		double x = patic.X[i] - L[0][0];
		double y = patic.Y[i] - L[1][0];
		double z = patic.Z[i] - L[2][0];
		int a = int(x*r_Grid[0] + eps);
		int b = int(y*r_Grid[1] + eps);
		int c = int(z*r_Grid[2] + eps);
		/*
		if ((a < 0) || (a >= B[0]) || (b < 0) || (b >= B[1]) || ((c < 0) || (c >= B[2])))
		{
			print("Patic %d outbound !!!!!!\n", i);
			continue;
		}
		*/
		if (a < 0) a = B[0] - 1;
		if (a >= B[0]) a = 0;
		if (b < 0) b = B[1] - 1;
		if (b >= B[1]) b = 0;
		if (c < 0) c = B[2] - 1;
		if (c >= B[2]) c = 0;

		int count = Locate3(a, b, c);
		bin[Locate4(a, b, c, binlen[count])] = i;
		binlen[count]++;
		binxid[i] = a;
		binyid[i] = b;
		binzid[i] = c;
	}
	printf("Updatebin done\n");
}

void UpdateNBList(Particle &patic)
{
	for (int i = 0; i<patic.Nr; i++)
	{
		int num = 0;
		for (int l = 0; l<VNum; l++)
		{
			int a = (binxid[i] + vect[l][0] + B[0]) % B[0];
			int b = (binyid[i] + vect[l][1] + B[1]) % B[1];
			int c = (binzid[i] + vect[l][2] + B[2]) % B[2];

			if (a < 0) a = B[0] - 1;
			if (a >= B[0]) a = 0;
			if (b < 0) b = B[1] - 1;
			if (b >= B[1]) b = 0;
			if (c < 0) c = B[2] - 1;
			if (c >= B[2]) c = 0;

			for (int k = 0; k<binlen[Locate3(a, b, c)]; k++)
			{
				int j = bin[Locate4(a, b, c, k)];
				double x = patic.X[i] - patic.X[j];
				double y = patic.Y[i] - patic.Y[j];
				double z = patic.Z[i] - patic.Z[j];
				double s = x*x + y*y + z*z;
				if ((i != j) && (s<dp_skin_2))
				{
					NBList[Locate2(num, i, NLen)] = j;
					NBDis[Locate2(num, i, NLen)] = dp - sqrt(s);
					num++;
				}
			}
		}
		NBLength[i] = num;
		if (num >= MinPts)
		{
			patic.Core[i] = 1;
			for (int h = 0; h < num; h++)
			{
				int number = NBList[Locate2(h, i, NLen)];
				patic.Border[number] = 1;
			}
		}
	}
	printf("Update NBL done.\n");
}

void assign_v2a(const vector<int>& clst, Cluster* iCluster, Particle &patic) // vector 2 array
{
	assert(iCluster);
	if (!clst.size())
		return;

	int Nr = clst.size();
	iCluster->size = Nr;
	for (int i = 0; i < Nr; i++) {
		iCluster->idList.push_back(clst[i]);
	}

fprintf(c_debug, "size=%d\n",iCluster->size);
fflush(c_debug);
	double cx = 0.0, cy = 0.0, cz = 0.0;

	// pre process
	int minx = 1, miny = 1, minz = 1;
	int maxx = -1, maxy = -1, maxz = -1;
	vector<int> xvec, yvec, zvec;
	vector<double> xpos, ypos, zpos;
	for (int i = 0; i < Nr; i++)
	{
		int j = clst[i];
		//gid_list[i] = patic_gid[j]; // update i_cluster->gid_list
		int xb = binxid[j];
		int yb = binyid[j];
		int zb = binzid[j];
		xvec.push_back(xb);
		yvec.push_back(yb);
		zvec.push_back(zb);
		xpos.push_back(patic.X[j]);
		ypos.push_back(patic.Y[j]);
		zpos.push_back(patic.Z[j]);
		minx = (minx > xb) ? xb : minx;
		miny = (miny > yb) ? yb : miny;
		minz = (minz > zb) ? zb : minz;
		maxx = (maxx > xb) ? maxx : xb;
		maxy = (maxy > yb) ? maxy : yb;
		maxz = (maxz > zb) ? maxz : zb;
	}

	int xcross = 0, ycross = 0, zcross = 0; // default no cross
	int findx, findy, findz;
	// x direction
	if (minx == 0 && maxx == B[0] - 1)
	{
		vector<int>::iterator it;
		for (findx = 1; findx < B[0] - 1; findx++)
		{
			it = find(xvec.begin(), xvec.end(), findx);
			if (it == xvec.end())
			{
				xcross = 1;
				break;
			}
		}
	}
	// y direction
	if (miny == 0 && maxy == B[1] - 1)
	{
		vector<int>::iterator it;
		for (findy = 1; findy < B[1] - 1; findy++)
		{
			it = find(yvec.begin(), yvec.end(), findy);
			if (it == yvec.end())
			{
				ycross = 1;
				break;
			}
		}
	}
	// z direction
	if (minz == 0 && maxz == B[2] - 1)
	{
		vector<int>::iterator it;
		for (findz = 1; findz < B[2] - 1; findz++)
		{
			it = find(zvec.begin(), zvec.end(), findz);
			if (it == zvec.end())
			{
				zcross = 1;
				break;
			}
		}
	}

	// update pos
	// x direction
	if (xcross)
	{
		for (int i = 0; i < Nr; i++)
		{
			if (xvec[i] < findx)
			{
				xpos[i] += PrsL[0];
			}
		}
	}
	// y direction
	if (ycross)
	{
		for (int i = 0; i < Nr; i++)
		{
			if (yvec[i] < findy)
			{
				ypos[i] += PrsL[1];
			}
		}
	}
	// z direction
	if (zcross)
	{
		for (int i = 0; i < Nr; i++)
		{
			if (zvec[i] < findz)
			{
				zpos[i] += PrsL[2];
			}
		}
	}

	// update centroid
	for (int i = 0; i < Nr; i++)
	{
		cx += xpos[i];
		cy += ypos[i];
		cz += zpos[i];
	}
	cx /= (Nr*1.0);
	cy /= (Nr*1.0);
	cz /= (Nr*1.0);

	// update the distance between every particle with centroid
	double tmprad2 = 0.0;
	for (int i = 0; i < Nr; i++)
	{
		double tmpx2 = (xpos[i] - cx)*(xpos[i] - cx);
		double tmpy2 = (ypos[i] - cy)*(ypos[i] - cy);
		double tmpz2 = (zpos[i] - cz)*(zpos[i] - cz);
		double tmpdis2 = tmpx2 + tmpy2 + tmpz2;
		tmprad2 += tmpdis2;
	}
	double tmprad = tmprad2 / (Nr*1.0);
	iCluster->diam = sqrt(tmprad) * 2;

	// debug
	if (iCluster->diam > 0.01 && Nr < 10)
	{
		fprintf(c_debug, "Nr=%d ", Nr);
		fflush(c_debug);
		for (int i = 0; i < Nr; i++)
		{
			int j = clst[i];
			fprintf(c_debug, "(%g %g %g) ", patic.X[j], patic.Y[j], patic.Z[j]);
			fflush(c_debug);
		}
		fprintf(c_debug, "\n");
		fflush(c_debug);
	}

	// update coordination number
	int sum = 0;
	for (int i = 0; i < Nr; i++)
	{
		int j = clst[i];
		sum += NBLength[j];
	}
	iCluster->CooNum = sum*1.0 / (Nr*1.0);

	printf("assign done.\n");
}

void DetectCluster(Particle &patic)
{
	int dc=0;
	int clstNum = 0;
	for (int i = 0; i < patic.Nr; i++)
	{
		if (!patic.Flag[i])
		{
			if (patic.Core[i])
			{
				patic.Flag[i] = 1;
				stack<int> point;
				vector<int> clst;
				point.push(i);
				clst.push_back(i);
			//	Cluster* iCluster = (Cluster*)malloc(sizeof(Cluster));
				Cluster* iCluster = new Cluster;
				while (!point.empty())
				{
					int k = point.top();
					point.pop();
					for (int j = 0; j < NBLength[k]; j++)
					{
						int tmp = NBList[Locate2(j, k, NLen)];
						if (tmp >= patic.Nr)
						{
							printf("Error!========================\n");
							exit(1);
						}
						if (patic.Flag[tmp]) continue;
						patic.Flag[tmp] = 1;
						clst.push_back(tmp);
						if (patic.Core[tmp])
							point.push(tmp);
					}
				}
fprintf(c_debug,"id=%d\t",clstNum);
fflush(c_debug);
				assign_v2a(clst, iCluster, patic);

				dbClst.push_back(iCluster);
				clstNum++;
printf("cluster num %d\n",clstNum);
				vector<int>(clst).swap(clst);
			}
			if (!patic.Core[i] && patic.Border[i])
				continue;
			if (!patic.Core[i] && !patic.Border[i])
				patic.Flag[i] = 1;
		}
	}
	ClusterNum = clstNum;
	printf("Detect Cluster done. Cluster num=%d\n",ClusterNum);
}

void OutputDBSCAN(Particle &patic, int filenum)
{
	FILE *fp;
	char dump[200] = "";
	char name[200];
	sprintf(name, "DBSCAN%d.vtk", filenum);
	fp = fopen(name, "w");

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET POLYDATA\n");
	fprintf(fp, "POINTS %d float\n", patic.Nr);

	for (int i = 0; i < patic.Nr; i++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\n", patic.X[i], patic.Y[i], patic.Z[i]);
	}
	fprintf(fp, "VERTICES %d %d\n", patic.Nr, patic.Nr * 2);
	for (int i = 0; i < patic.Nr; i++)
		fprintf(fp, "%5d\t%8d\n", 1, i);

	fprintf(fp, "POINT_DATA %d\n", patic.Nr);
	fprintf(fp, "SCALARS clusterid float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	// update particle's inClst
	for (int i = 0; i < ClusterNum; i++) {
		Cluster* ic = dbClst[i];
		for (int j = 0; j < ic->size; j++) {
			int pid = ic->idList[j];
			patic.inClst[pid] = i + 1;
		}
	}
	// output clusterid
	for (int i = 0; i < patic.Nr; i++)
		fprintf(fp, "%d\n", patic.inClst[i]);

	fclose(fp);
	printf("Output DBSCAN cluster done.\n");
}

void OutputResult(int filenum)
{
	char name[200];
	sprintf(name,"dbresult-%.3f-%d-%.3f.csv",FileTime,MinPts,Skin);
	FILE* fp=fopen(name,"w");
	fprintf(fp,"clusterid,np,diam,coonum\n");
	for(int i=0;i<dbClst.size();i++) {
		Cluster* ic=dbClst[i];
		fprintf(fp,"%d,%d,%2.6e,%2.6e\n",i,ic->size,ic->diam,ic->CooNum);
	}
	fclose(fp);
	printf("Output result done\n");
}

void OneStep()
{
	ReadDBSCANInfo();
	InitializeNeighbourList();
	char filename[200];
	int filenum = int(FileTime*scale + eps);
	sprintf(filename, "part%d.vtk", filenum);
	ReadVTK(filename, p);
	Updatebin(0, p.Nr, p);
	UpdateNBList(p);
	DetectCluster(p);
	if (outputVTK) {
		OutputDBSCAN(p,filenum);
	}
	OutputResult(filenum);
}

int main()
{
	c_debug = fopen("ClusterDebug.dat", "w");
	OneStep();
	fclose(c_debug);

	return 0;
}

