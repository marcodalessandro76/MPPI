"""
This module implement the merge of the content of a list of ndb.QP databases
into a single database.
The function is deeply inspired from the merge_qp function of yambopy
"""

import numpy as np
from netCDF4 import Dataset


# def parse_ndbQP(file):
#     """
#     Parse the content of the ndb.QP
#
#     Args:
#         file (:py:class:`string`) : name of the file with the database, including the path
#
#     Return:
#         :py:class:`list` :.....
#     """
#
#     dataset = Dataset(file)
#     try:
#         qp_test = datasets[0]['QP_E']
#     except IndexError:
#         print('Old version of database detected')
#     PARS = list(map(int,d['PARS'][:]))
#     nkpoints, nqps, nstrings = PARS[1],PARS[2],PARS[-1]
#     QP_table = d['QP_table'][:].T
#     QP_kpts = d['QP_kpts'][:].T
#     QP_E = d['QP_E'][:]
#     QP_E0 = d['QP_Eo'][:]
#     QP_Z = d['QP_Z'][:]

def merge_qp(output,files):
    """
    Merge the quasiparticle databases produced by yambo
    
    Args:
        output (:py:class:`string`) : name of the output file
        files  (:py:class:`list`)   : list of the input files to be merged

    """
    #read all the files and display main info in each of them
    #filenames = [ f.name for f in files]
    #datasets  = [ Dataset(filename) for filename in filenames]
    datasets  = [ Dataset(f) for f in files]
    
    #call compatibility version if old dataset detected
    try:
        qp_test = datasets[0]['QP_E']
    except IndexError:
        print('Old version of database detected')
        raise IndexError('Problem with the databases')
    #finally:

    print("=========input=========")
    QP_table, QP_kpts, QP_E, QP_E0, QP_Z = [], [], [], [], []
    #for d,filename in zip(datasets,filenames):
    for d,filename in zip(datasets,files):    
        #PARS = list(map(int,d['PARS'][:]))
        PARS = list(map(int, np.ma.filled(d['PARS'][:], 0)))
        nkpoints, nqps, nstrings = PARS[1],PARS[2],PARS[-1]
        #_, nkpoints, nqps, _, nstrings = list(map(int,d['PARS'][:]))
        print("filename:    ", filename)
        QP_table.append( d['QP_table'][:].T )
        QP_kpts.append( d['QP_kpts'][:].T )
        QP_E.append( d['QP_E'][:] )
        QP_E0.append( d['QP_Eo'][:] )
        QP_Z.append( d['QP_Z'][:] )

    # create the QP_table
    QP_table_save = np.vstack(QP_table)

    # create the kpoints table
    #create a list with the bigger size of QP_table
    nkpoints = int(max(QP_table_save[:,2]))
    QP_kpts_save = np.zeros([nkpoints,3])
    #iterate over the QP's and store the corresponding kpoint
    for qp_file,kpts in zip(QP_table,QP_kpts):
        #iterate over the kpoints and save the coordinates on the list
        for qp in qp_file:
            try:               n1,n2,nk = list(map(int,qp))
            except ValueError: n1,n2,nk,ns = list(map(int,qp))
            QP_kpts_save[nk-1] = kpts[nk-1]

    # create the QPs energies table
    QP_E_save  = np.concatenate(QP_E,axis=0)
    QP_E0_save = np.concatenate(QP_E0)
    QP_Z_save  = np.concatenate(QP_Z,axis=0)

    #create reference file from one of the files
    netcdf_format = datasets[0].data_model
    fin  = datasets[0]
    fout = Dataset(output,'w',format=netcdf_format)

    variables_update = ['QP_table', 'QP_kpts', 'QP_E', 'QP_Eo', 'QP_Z']
    variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_save, QP_E0_save, QP_Z_save]
    variables_dict   = dict(list(zip(variables_update,variables_save)))
    PARS_save = fin['PARS'][:]
    PARS_save[1:3] = nkpoints,len(QP_table_save)

    #create the description string
    kmin,kmax = np.amin(QP_table_save[:,2]),np.amax(QP_table_save[:,2])
    bmin,bmax = np.amin(QP_table_save[:,1]),np.amax(QP_table_save[:,1])
    description = "QP @ K %03d - %03d : b %03d - %03d"%(kmin,kmax,bmin,bmax)
    description_save = np.array([i for i in " %s"%description])
    QP_k_range, QP_b_range = [kmin,kmax], [bmin,bmax]


    #output data
    print("========output=========")
    print("filename:    ", output)
    print("description: ", description)

    #copy dimensions
    for dname, the_dim in list(fin.dimensions.items()):
        fout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    #get dimensions
    def dimensions(array):
        return tuple([ 'D_%010d'%d for d in array.shape ])

    #create missing dimensions
    for v in variables_save:
        for dname,d in zip( dimensions(v),v.shape ):
            if dname not in list(fout.dimensions.keys()):
                fout.createDimension(dname, d)

    #copy variables
    for v_name, varin in list(fin.variables.items()):
        if v_name in variables_update:
            #get the variable
            merged = variables_dict[v_name]
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, dimensions(merged))
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            #save outvar
            outVar[:] = merged

        else:
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, varin.dimensions)
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            if v_name=='PARS':
                outVar[:] = PARS_save[:]
            elif v_name=='DESC_strings_%05d'%(nstrings):
                outVar[:] = varin[:]
                outVar[:,:len(description_save)] = description_save.T
            elif v_name=='QP_QP_@_state_1_K_range':
                outVar[:]=QP_k_range
            elif v_name=='QP_QP_@_state_1_b_range':
                outVar[:]=QP_b_range
            else:
                outVar[:] = varin[:]

    fout.close()
