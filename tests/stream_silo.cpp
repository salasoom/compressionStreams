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

template <typename T> static void StreamCoordsOut(std::ostream &ostr,
    int ndims, int npts, T const * const *coords)
{
    for (int i = 0; i < npts; i++)
    {
        ostr << coords[0][i];
        for (int j = 1; j < ndims; j++)
            ostr << " " << coords[j][i];
        ostr << std::endl;
    }
}

template <typename T> static void StreamCoordsIn(std::stringstream &istr,
    int ndims, int npts, T **coords)
{
    for (int j = 0; j < ndims; j++)
        coords[j] = (T*) malloc(npts * sizeof(T));

    for (int i = 0; i < npts; i++)
    {
        for (int j = 0; j < ndims; j++)
            istr >> coords[j][i];
    }
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

template <typename T> static bool StreamKeyValueIn(std::stringstream &istrstrm,
    size_t &n, std::string const &key, T& val)
{
    std::string const &str = istrstrm.str();
    n = shortfind(str, key, n);
    if (n == str.npos)
        return false;
    istrstrm.seekg(n+key.size());
    istrstrm >> val;
    return true;
}

static void StreamPointmeshOut(DBpointmesh const *m, std::ostream &ostr)
{
    // output floating point coordinates
    ostr << "nnodes = " << m->nels << std::endl;
    ostr << "ndims = " << m->ndims << std::endl;
    ostr << "datatype = " << (m->datatype == DB_DOUBLE ? "double" : "float") << std::endl;
    if (m->datatype == DB_DOUBLE)
    {
        ostr << std::setprecision(16);
        StreamCoordsOut(ostr, m->ndims, m->nels, (double const * const *)m->coords);
    }
    else
    {
        ostr << std::setprecision(8);
        StreamCoordsOut(ostr, m->ndims, m->nels, (float const * const *)m->coords);
    }
}

static void StreamQuadmeshOut(DBquadmesh const *m, std::ostream &ostr)
{
    ostr << "ndims = " << m->ndims << std::endl;
    ostr << "nnodes = " << m->nnodes << std::endl;
    ostr << "dims = " << m->dims[0];
    for (int j = 1; j < m->ndims; j++)
        ostr << " " << m->dims[j];
    ostr << std::endl;
    ostr << "datatype = " << (m->datatype == DB_DOUBLE ? "double" : "float") << std::endl;
    if (m->coordtype == DB_COLLINEAR)
    {
        ostr << "coordtype = collinear" << std::endl;
        if (m->datatype == DB_DOUBLE)
        {
            ostr << std::setprecision(16);
            for (int j = 0; j < m->ndims; j++)
            {
                double const *coords = (double const *) m->coords[j];
                StreamCoordsOut(ostr, 1, m->dims[j], &coords);
            }
        }
        else
        {
            ostr << std::setprecision(8);
            for (int j = 0; j < m->ndims; j++)
            {
                float const *coords = (float const *) m->coords[j];
                StreamCoordsOut(ostr, 1, m->dims[j], &coords);
            }
        }
    }
    else
    {
        ostr << "coordtype = non-collinear" << std::endl;
        int npts = 1;
        for (int j = 0; j < m->ndims; j++, npts *= m->dims[j]);
        if (m->datatype == DB_DOUBLE)
        {
            ostr << std::setprecision(16);
            StreamCoordsOut(ostr, m->ndims, npts, (double const * const *) m->coords);
        }
        else
        {
            ostr << std::setprecision(8);
            StreamCoordsOut(ostr, m->ndims, npts, (float const * const *) m->coords);
        }
    }
}

static void StreamUcdmeshOut(DBucdmesh const *m, std::ostream &ostr)
{
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
        ostr << std::setprecision(16);
        StreamCoordsOut(ostr, m->ndims, m->nnodes, (double const * const *)m->coords);
    }
    else
    {
        ostr << std::setprecision(8);
        StreamCoordsOut(ostr, m->ndims, m->nnodes, (float const * const *)m->coords);
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

static DBpointmesh *StreamPointmeshIn(std::stringstream &istrstrm)
{
    std::string const &str = istrstrm.str();
    DBpointmesh *pm = DBAllocPointmesh();
    std::string dtstr;
    size_t n = 0;

    if (!StreamKeyValueIn(istrstrm, n, "nnodes = ", pm->nels)) return 0;
    if (!StreamKeyValueIn(istrstrm, n, "ndims = ", pm->ndims)) return 0;
    if (!StreamKeyValueIn(istrstrm, n, "datatype = ", dtstr)) return 0;
    pm->datatype = (dtstr == "double" ? DB_DOUBLE : DB_FLOAT);

    if (pm->datatype == DB_DOUBLE)
        StreamCoordsIn(istrstrm, pm->ndims, pm->nels, (double **) pm->coords);
    else
        StreamCoordsIn(istrstrm, pm->ndims, pm->nels, (float **) pm->coords);

    return pm;
}

static DBquadmesh *StreamQuadmeshIn(std::stringstream &istrstrm)
{
    std::string const &str = istrstrm.str();
    DBquadmesh *qm = DBAllocQuadmesh();
    std::string tstr;
    size_t n = 0;

    if (!StreamKeyValueIn(istrstrm, n, "ndims = ", qm->ndims)) return 0;
    if (!StreamKeyValueIn(istrstrm, n, "nnodes = ", qm->nnodes)) return 0;
    if (!StreamKeyValueIn(istrstrm, n, "dims = ", qm->dims[0])) return 0;
    for (int j = 1; j < qm->ndims; j++)
        istrstrm >> qm->dims[j];
    if (!StreamKeyValueIn(istrstrm, n, "datatype = ", tstr)) return 0;
    qm->datatype = (tstr == "double" ? DB_DOUBLE : DB_FLOAT);
    if (!StreamKeyValueIn(istrstrm, n, "coordtype = ", tstr)) return 0;
    qm->coordtype = (tstr == "collinear" ? DB_COLLINEAR : DB_NONCOLLINEAR);
    if (qm->coordtype == DB_COLLINEAR)
    {
        if (qm->datatype == DB_DOUBLE)
        {
            for (int j = 0; j < qm->ndims; j++)
                StreamCoordsIn(istrstrm, 1, qm->dims[j], (double**) &(qm->coords[j]));
        }
        else
        {
            for (int j = 0; j < qm->ndims; j++)
                StreamCoordsIn(istrstrm, 1, qm->dims[j], (float**) &(qm->coords[j]));
        }
    }
    else
    {
        if (qm->datatype == DB_DOUBLE)
            StreamCoordsIn(istrstrm, qm->ndims, qm->nnodes, (double **) qm->coords);
        else
            StreamCoordsIn(istrstrm, qm->ndims, qm->nnodes, (float **) qm->coords);
    }

    return qm;
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

template <typename T> static bool CompareArrays(int ndims, int npts,
    T const * const *coordsa, T const * const *coordsb)
{
    for (int j = 0; j < ndims; j++)
    {
        for (int i = 0; i < npts; i++)
        {
            if (coordsa[j][i] != coordsb[j][i])
                return false;
        }
    }
    return true;
}

static bool ComparePointmesh(DBpointmesh const *pma, DBpointmesh *pmb)
{
    if (pma->nels != pmb->nels) return false;
    if (pma->ndims != pmb->ndims) return false;
    if (pma->datatype != pmb->datatype) return false;
    if (pma->datatype == DB_DOUBLE)
    {
        if (!CompareArrays(pma->ndims, pma->nels,
             (double const * const *) pma->coords, (double const * const *) pmb->coords))
            return false;
    }
    else
    {
        if (!CompareArrays(pma->ndims, pma->nels,
             (float const * const *) pma->coords, (float const * const *) pmb->coords))
            return false;
    }

    return true;
}

static bool CompareQuadmesh(DBquadmesh const *qma, DBquadmesh *qmb)
{
    if (qma->nnodes != qmb->nnodes) return false;
    if (qma->ndims != qmb->ndims) return false;
#if 0
    if (!CompareArrays(1, qma->ndims,
         (int const * const *) &(qma->dims[0]),
         (int const * const *) &(qmb->dims[0])))
        return false;
#endif
    if (qma->datatype != qmb->datatype) return false;
    if (qma->coordtype != qmb->coordtype) return false;
    if (qma->coordtype == DB_COLLINEAR)
    {
        if (qma->datatype == DB_DOUBLE)
        {
            for (int j = 0; j < qma->ndims; j++)
            {
                if (!CompareArrays(1, qma->dims[j],
                     (double const * const *)&(qma->coords[j]),
                     (double const * const *)&(qmb->coords[j])))
                    return false;
            } 
        }
        else
        {
            for (int j = 0; j < qma->ndims; j++)
            {
                if (!CompareArrays(1, qma->dims[j],
                     (float const * const *)&(qma->coords[j]),
                     (float const * const *)&(qmb->coords[j])))
                    return false;
            } 
        }
    }
    else
    {
        if (qma->datatype == DB_DOUBLE)
        {
            if (!CompareArrays(qma->ndims, qma->nnodes,
                 (double const * const *)qma->coords,
                 (double const * const *)qmb->coords))
                return false;
        }
        else
        {
            if (!CompareArrays(qma->ndims, qma->nnodes,
                 (float const * const *)qma->coords,
                 (float const * const *)qmb->coords))
                return false;
        }
    }

    return true;
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
