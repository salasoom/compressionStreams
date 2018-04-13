#include <silo.h>

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iomanip>

#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>
#include <string>

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
    // output floating point coordinates
    ostr << "nnodes = " << m->nels << std::endl;
    ostr << "ndims = " << m->ndims << std::endl;
    ostr << "datatype = " << (m->datatype == DB_DOUBLE ? "double" : "float") << std::endl;
    if (m->datatype == DB_DOUBLE)
    {
        double *coords[3] = {(double *) m->coords[0],
                             (double *) m->coords[1],
                             (double *) m->coords[2]};
        ostr << std::setprecision(16);
        for (int i = 0; i < m->nels; i++)
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
        ostr << std::setprecision(8);
        for (int i = 0; i < m->nels; i++)
        {
            ostr << coords[0][i];
            for (int j = 1; j < m->ndims; j++)
                ostr << " " << coords[j][i];
            ostr << std::endl;
        }
    }
}

static void StreamQuadmeshOut(DBquadmesh const *m, std::ostream &ostr)
{
    ostr << "nnodes = " << m->nnodes << std::endl;
    ostr << "ndims = " << m->ndims << std::endl;
    ostr << "dims = " << m->dims[0];
    ostr << "datatype = " << (m->datatype == DB_DOUBLE ? "double" : "float") << std::endl;
    for (int i = 1; i < m->ndims; i++)
        ostr << " " << m->dims[i];
    if (m->coordtype == DB_COLLINEAR)
    {
        ostr << "coordtype = collinear" << std::endl;
        if (m->datatype == DB_DOUBLE)
        {
        }
        else
        {
        }
    }
    else
    {
        ostr << "coordtype = non-collinear" << std::endl;
        if (m->datatype == DB_DOUBLE)
        {
        }
        else
        {
        }
    }
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
        ostr << "nzones = " << zl->nzones << std::endl;
        ostr << "nshapes = " << zl->nshapes << std::endl;
        for (int i = 0; i < zl->nshapes; i++)
        {
            ostr << "shapesize = " << zl->shapesize[i] << std::endl;
            ostr << "shapetype = " << zl->shapetype[i] << std::endl;
            ostr << "shapecount = " << zl->shapecnt[i] << std::endl;
            for (int j = 0; j < zl->shapecnt[i]; j++)
            {
                ostr << zl->nodelist[n++];
                for (int k = 1; k < zl->shapesize[i]; k++)
                    ostr << " " << zl->nodelist[n++];
                ostr << std::endl;
            }
        }
    }

    // output floating point coordinates
    ostr << "nnodes = " << m->nnodes << std::endl;
    ostr << "ndims = " << m->ndims << std::endl;
    ostr << "datatype = " << (m->datatype == DB_DOUBLE ? "double" : "float") << std::endl;
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

// Like a normal string::find() but limit search to maxlen chars
static size_t shortfind(std::string const &str, std::string const &needle,
    size_t start, size_t maxlen=100)
{
    std::string const &shortstr = str.substr(start, maxlen);
    size_t n = shortstr.find(needle);
    if (n == shortstr.npos)
        return str.npos;
    return start + n;
}

static DBpointmesh *StreamPointmeshIn(std::stringstream &istrstrm)
{
    std::string const &str = istrstrm.str();
    DBpointmesh *pm = DBAllocPointmesh();
    size_t n;

    // read nnodes
    std::string tagstr = "nnodes = ";
    n = shortfind(str, tagstr, 0);
    if (n == str.npos)
        return 0;
    istrstrm.seekg(n+tagstr.size());
    istrstrm >> pm->nels;

    tagstr = "ndims = ";
    n = shortfind(str, tagstr, n);
    if (n == str.npos)
        return 0;
    istrstrm.seekg(n+tagstr.size());
    istrstrm >> pm->ndims;

    tagstr = "datatype = ";
    n = shortfind(str, tagstr, n);
    if (n == str.npos)
        return 0;
    istrstrm.seekg(n+tagstr.size());
    istrstrm >> tagstr;
    pm->datatype = (tagstr == "double" ? DB_DOUBLE : DB_FLOAT);

    if (pm->datatype == DB_DOUBLE)
    {
        for (int i = 0; i < pm->ndims; i++)
            pm->coords[i] = (double *) malloc(pm->nels * sizeof(double));
        for (int i = pm->ndims; i < 3; i++)
            pm->coords[i] = 0;
        double *coords[3] = {(double*)pm->coords[0],
                             (double*)pm->coords[1],
                             (double*)pm->coords[2]};
        for (int i = 0; i < pm->nels; i++)
        {
            for (int j = 0; j < pm->ndims; j++)
                istrstrm >> coords[j][i];
        }
    }
    else
    {
        for (int i = 0; i < pm->ndims; i++)
            pm->coords[i] = (float *) malloc(pm->nels * sizeof(float));
        for (int i = pm->ndims; i < 3; i++)
            pm->coords[i] = 0;
        float *coords[3] = {(float*)pm->coords[0],
                            (float*)pm->coords[1],
                            (float*)pm->coords[2]};
        for (int i = 0; i < pm->nels; i++)
        {
            for (int j = 0; j < pm->ndims; j++)
                istrstrm >> coords[j][i];
        }
    }

    return pm;
}

static DBquadmesh *StreamQuadmeshIn(std::stringstream &istrstrm)
{
    return 0;
}

static DBucdmesh *StreamUcdmeshIn(std::stringstream &istrstrm)
{
    return 0;
}

static void *StreamMeshIn(int mt, char const *ifile)
{
    // open the input stream
    std::ifstream meshfile(ifile);
    std::stringstream meshfilestr;
    meshfilestr << meshfile.rdbuf(); // not good if large files

    if (mt == DB_POINTMESH)
        return StreamPointmeshIn(meshfilestr);
    else if (mt == DB_UCDMESH)
        return StreamUcdmeshIn(meshfilestr);
    else
        return StreamQuadmeshIn(meshfilestr);
}

static bool ComparePointmesh(DBpointmesh const *pma, DBpointmesh *pmb)
{
    if (pma->nels != pmb->nels) return false;
    if (pma->ndims != pmb->ndims) return false;
    if (pma->datatype != pmb->datatype) return false;
    if (pma->datatype == DB_DOUBLE)
    {
        double const *coordsa[3] = {(double*)pma->coords[0],
                                    (double*)pma->coords[1],
                                    (double*)pma->coords[2]};
        double const *coordsb[3] = {(double*)pmb->coords[0],
                                    (double*)pmb->coords[1],
                                    (double*)pmb->coords[2]};
        for (int j = 0; j < pma->ndims; j++)
        {
            for (int i = 0; i < pma->nels; i++)
            {
                if (coordsa[j][i] != coordsb[j][i]) return false;
            }
        }
    }
    else
    {
        float const *coordsa[3] = {(float*)pma->coords[0],
                                   (float*)pma->coords[1],
                                   (float*)pma->coords[2]};
        float const *coordsb[3] = {(float*)pmb->coords[0],
                                   (float*)pmb->coords[1],
                                   (float*)pmb->coords[2]};
        for (int j = 0; j < pma->ndims; j++)
        {
            for (int i = 0; i < pma->nels; i++)
            {
                if (coordsa[j][i] != coordsb[j][i]) return false;
            }
        }
    }

    return true;
}

static bool CompareQuadmesh(DBquadmesh const *qma, DBquadmesh *qmb)
{
    return false;
}

static bool CompareUcdmesh(DBucdmesh const *uma, DBucdmesh *umb)
{
    return false;
}

static bool CompareMesh(int mt, void const *mesh, void const *zmesh)
{
    if (mt == DB_POINTMESH)
        return ComparePointmesh((DBpointmesh*)mesh,(DBpointmesh*)zmesh);
    else if (mt == DB_UCDMESH)
        return CompareUcdmesh((DBucdmesh*)mesh,(DBucdmesh*)zmesh);
    else
        return CompareQuadmesh((DBquadmesh*)mesh,(DBquadmesh*)zmesh);
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

    strcpy(ofile, "stream_silo_out.txt");
    strcpy(ifile, "../data/globe_hdf5.silo");
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

    // Stream the mesh back in
    void *zmesh = StreamMeshIn(mt, ofile);

    /* compare streamd out/in mesh to original */
    if (CompareMesh(mt, mesh, zmesh))
    {
        std::cout << "Meshes are identical" << std::endl;
    }
    else
    {
        std::cout << "Meshes FAILED" << std::endl;
    }

    return 0;
}
