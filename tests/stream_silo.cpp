#include <silo.h>

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iomanip>

#include <fstream>

#define NAME_LEN 256

/* convenience macro to handle command-line args and help */
#define HANDLE_ARG(A,PARSEA,PRINTA,HELPSTR)                     \
{                                                               \
    int i;                                                      \
    char tmpstr[64];                                            \
    int len;                                                    \
    int len2 = strlen(#A)+1;                                    \
    for (i = 0; i < argc; i++)                                  \
    {                                                           \
        if (!strncmp(argv[i], #A"=", len2))                     \
        {                                                       \
            A = PARSEA;                                         \
            if (!strncasecmp(argv[i], "help", 4))               \
            {                                                   \
                fflush(stdout);                                 \
                return 0;                                       \
            }                                                   \
            else                                                \
            {                                                   \
                break;                                          \
            }                                                   \
        }                                                       \
    }                                                           \
    len = snprintf(tmpstr, sizeof(tmpstr), "%s=" PRINTA, #A, A);\
    printf("    %s%*s\n",tmpstr,60-len,#HELPSTR);               \
}


/* convenience macro to handle errors */
#define ERROR(FNAME)                                              \
do {                                                              \
    int _errno = errno;                                           \
    fprintf(stderr, #FNAME " failed at line %d, errno=%d (%s)\n", \
        __LINE__, _errno, _errno?strerror(_errno):"ok");          \
    return 1;                                                     \
} while(0)

static void *
GetMesh(DBfile *dbfile, char const *mname, int &mt)
{
    mt = DBInqMeshtype(dbfile, mname);

    // Get the actual mesh object
    if (mt == DB_POINTMESH)
        return DBGetPointmesh(dbfile, mname);
    else if (mt == DB_UCDMESH)
        return DBGetUcdmesh(dbfile, mname);
    else
        return DBGetQuadmesh(dbfile, mname);
}

static void StreamPointmeshOut(DBpointmesh const *m, std::ostream &ostr)
{
}
static void StreamQuadmeshOut(DBquadmesh const *m, std::ostream &ostr)
{
}
static void StreamUcdmeshOut(DBucdmesh const *m, std::ostream &ostr)
{

    int           *shapecnt;    /* [nshapes] occurences of each shape */
    int           *shapesize;   /* [nshapes] Number of nodes per shape */
    int           *shapetype;   /* [nshapes] Type of shape */
    int           *nodelist;    /* Sequent lst of nodes which comprise zones */
    int            lnodelist;   /* Number of nodes in nodelist */
    int            origin;      /* '0' or '1' */

    // output integer topology
    if (m->zones)
    {
        int n = 0;
        DBzonelist *zl = m->zones;
        for (int i = 0; i < zl->nshapes; i++)
        {
            for (int j = 0; j < zl->shapecnt[i]; j++)
            {
                ostr << zl->shapesize[i];
                for (int k = 0; k < zl->shapesize[i]; k++)
                    ostr << " " << zl->nodelist[n++];
                ostr << std::endl;
            }
        }
    }

    // output floating point coordinates
    if (m->datatype == DB_DOUBLE)
    {
        double *coords[3] = {(double *) m->coords[0],
                             (double *) m->coords[1],
                             (double *) m->coords[2]};
        ostr << std::setprecision(14);
        for (int i = 0; i < m->nnodes; i++)
        {
            ostr << coords[0][i];
            for (int j = 1; j < m->ndims; j++)
                ostr << " " << coords[j][i];
            ostr << std::endl;
        }
    }
    else
    {
        float *coords[3] = {(float *) m->coords[0],
                            (float *) m->coords[1],
                            (float *) m->coords[2]};
        ostr << std::setprecision(7);
        for (int i = 0; i < m->nnodes; i++)
        {
            ostr << coords[0][i];
            for (int j = 1; j < m->ndims; j++)
                ostr << " " << coords[j][i];
            ostr << std::endl;
        }
    }
}

static void
StreamMeshOut(int mt, void *mesh, char const *ofile, int gzlevel)
{
    // open the outstream
    std::ofstream meshfile(ofile);

    if (mt == DB_POINTMESH)
        StreamPointmeshOut((DBpointmesh*)mesh, meshfile);
    else if (mt == DB_UCDMESH)
        StreamUcdmeshOut((DBucdmesh*)mesh, meshfile);
    else
        StreamQuadmeshOut((DBquadmesh*)mesh, meshfile);
}

int main(int argc, char **argv)
{

    int i;

    char *ifile = (char *) calloc(NAME_LEN,sizeof(char));
    char *mname = (char *) calloc(NAME_LEN,sizeof(char));
    char *ofile = (char *) calloc(NAME_LEN,sizeof(char));
    char *zalg = (char *) calloc(NAME_LEN,sizeof(char));

    /* compression parameters (defaults taken from ZFP header) */
    int help = 0;
    int gzlevel = 5;
    int zfpmode = 3;
    double rate = 4;
    double acc = 0;
    uint prec = 11;
    uint minbits = 0;
    uint maxbits = 4171;
    uint maxprec = 64;
    int minexp = -1074;
    int *ibuf = 0;
    double *buf = 0;

    HANDLE_ARG(ifile,strndup(argv[i]+len2,NAME_LEN), "\"%s\"",set input Silo filename);
    HANDLE_ARG(mname,strndup(argv[i]+len2,NAME_LEN), "\"%s\"",set input Mesh name);
    HANDLE_ARG(ofile,strndup(argv[i]+len2,NAME_LEN), "\"%s\"",set output stream filename);
    HANDLE_ARG(zalg,strndup(argv[i]+len2,NAME_LEN), "\"%s\"",set compression algorithm);
    HANDLE_ARG(gzlevel,(int) strtol(argv[i]+len2,0,10),"%d",set gzip level (1-9));
    HANDLE_ARG(zfpmode,(int) strtol(argv[i]+len2,0,10),"%d",set zfp mode (1=rate,2=prec,3=acc,4=expert));
    HANDLE_ARG(rate,(double) strtod(argv[i]+len2,0),"%g",set rate for rate mode of filter);
    HANDLE_ARG(acc,(double) strtod(argv[i]+len2,0),"%g",set accuracy for accuracy mode of filter);
    HANDLE_ARG(prec,(uint) strtol(argv[i]+len2,0,10),"%u",set precision for precision mode of zfp filter);
    HANDLE_ARG(minbits,(uint) strtol(argv[i]+len2,0,10),"%u",set minbits for expert mode of zfp filter);
    HANDLE_ARG(maxbits,(uint) strtol(argv[i]+len2,0,10),"%u",set maxbits for expert mode of zfp filter);
    HANDLE_ARG(maxprec,(uint) strtol(argv[i]+len2,0,10),"%u",set maxprec for expert mode of zfp filter);
    HANDLE_ARG(minexp,(int) strtol(argv[i]+len2,0,10),"%d",set minexp for expert mode of zfp filter);
    HANDLE_ARG(help,(int)strtol(argv[i]+len2,0,10),"%d",this help message);

    // Open the Silo file
    DBfile *dbfile = DBOpen(ifile, DB_UNKNOWN, DB_READ); 

    int mt; void *mesh = GetMesh(dbfile, mname, mt);

    // Close the silo file
    DBClose(dbfile);

    // Stream out the mesh
    // just level for now
    StreamMeshOut(mt, mesh, ofile, gzlevel);

#if 0
    // Stream the mesh back in
    StreamMeshIn(mt, zmesh, ofile);

    /* compare streamd out/in mesh to original */
    CompareMesh(mt, mesh, zmesh);
#endif

    return 0;
}
