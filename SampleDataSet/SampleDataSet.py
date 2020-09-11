from paraview.util.vtkAlgorithm import *

#------------------------------------------------------------------------------
# A filter example.
#------------------------------------------------------------------------------
@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)

class SampleDataSet(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.vars = ['aaa', 'bbb']
        self.minSpacing = 1.0
        self.maxSpacing = 0.1
        self.samples = 1000
        self.retained = None

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
        return 1

    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return self.vars

    @smproperty.stringvector(name="String", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="StringInfo" function="StringInfo"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetString(self, value):
        self.inputArray = value
        self.Modified()
        return 

    @smproperty.intvector(name="Continuity", label="Continuity", default_values=1)
    @smdomain.intrange(min=0, max=1)
    def SetContinuity(self, x):
        self.continuity = x
        self.Modified()

    @smproperty.intvector(name="Samples", label="Samples", default_values=1000)
    @smdomain.intrange(min=0, max=1)
    def SetSamples(self, x):
        self.samples = x
        self.Modified()

    @smproperty.doublevector(name="minSpacing", label="spacing at min probability", default_values=4.0)
    @smdomain.doublerange(min=0, max=1000000)
    def SetMinSpacing(self, x):
        self.minSpacing = x
        self.Modified()

    @smproperty.doublevector(name="maxSpacing", label="spacing at max probability", default_values=1.0)
    @smdomain.doublerange(min=0, max=1000000)
    def SetMaxSpacing(self, x):
        self.maxSpacing = x
        self.Modified()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkDataSet
        executive = self.GetExecutive()
        input0 = vtkDataSet.GetData(inInfoVec[0], 0)

        vars = []
        for i in range(input0.GetPointData().GetNumberOfArrays()):
          vars.append(input0.GetPointData().GetArray(i).GetName())

        if 'pdf' in vars:
          self.vars = ['pdf'] 
        else:
          self.vars = []
        for i in vars:
          if i != 'pdf':
            self.vars.append(i)
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        import sys, os
        import vtk
        for i in dir(vtk):
          if i[:5] == 'vtkRe':
            print(i)
        from vtk.numpy_interface import dataset_adapter as dsa
        import ctypes as C
        import numpy as np

        if 'HOME' in os.environ:
          ccode_so = os.environ['HOME'] + '/ParaviewPythonModules/libSampleDataset.so'
        elif 'HOMEPATH' in os.environ:
          ccode_so = os.environ['HOMEPATH'] + '/ParaviewPythonModules/libSampleDataset.so'
        else:
          ccode_so =  'SampleDataset.so'

        if not os.path.isfile(ccode_so):
          print('can\'t find so:', ccode_so)
          return

        ccode = C.CDLL(ccode_so)
        if not ccode:
          print('failed to load ', ccode.so)

        ccode.SampleCartesian.argtypes = [C.c_int,                                                    # desired number of samples
                                          np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # origin  (3)
                                          np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # spacing (3)
                                          np.ctypeslib.ndpointer(C.c_int, flags="C_CONTIGUOUS"),      # counts  (3)
                                          C.c_float,                                                  # minimum spacing
                                          C.c_float,                                                  # maximum spacing
                                          C.c_int,                                                    # pdep data?
                                          np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # data
                                          C.c_int,                                                    # number of retained samples
                                          np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # retained samples
                                          np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS")]    # new data interpolated on retained samples

        ccode.SampleTetrahedra.argtypes = [C.c_int,                                                    # desired number of samples
                                           C.c_int,                                                    # number of cells
                                           C.c_float,                                                  # minimum spacing
                                           C.c_float,                                                  # maximum spacing
                                           np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # points
                                           np.ctypeslib.ndpointer(C.c_int, flags="C_CONTIGUOUS"),      # cells
                                           C.c_int,                                                    # pdep data?
                                           np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # data
                                           C.c_int,                                                    # number of retained samples
                                           np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),    # retained samples
                                           np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS")]    # new data interpolated on retained samples

        ccode.GetNumberOfSamples.restype = C.c_int;
        ccode.GetSamples.argtypes = [np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS")];

        # self.inputArray = 'PDF'
        print("Sample USING ", self.inputArray)

        if self.retained != None:
          r = vtk.vtkResampleWithDataSet()
          r.SetSourceData(vtk.vtkDataSet.GetData(inInfoVec[0], 0))
          r.SetInputData(self.retained)
          r.Update()
          rr = dsa.WrapDataObject(r.GetOutput())
          del r
          rs = rr.Points.astype('f4')
          rd = rr.PointData[self.inputArray].astype('f4')
          nRetained = rr.GetNumberOfPoints()
          del rr
        else:
          nRetained = 0
          rs = np.zeros(1).astype('f4')
          rd = np.zeros(1).astype('f4')

        inpt = vtk.vtkImageData.GetData(inInfoVec[0], 0)
        if inpt != None:
          inpt = dsa.WrapDataObject(inpt)

          o = inpt.VTKObject.GetOrigin()
          e = inpt.VTKObject.GetExtent()
          s = inpt.VTKObject.GetSpacing()

          k = np.array([e[i*2+1] - e[i*2] + 1 for i in range(3)]).astype('i4')
          o = np.array([o[i] + e[2*i]*s[i] for i in range(3)]).astype('f4')
          s = np.array(s).astype('f4')

          data    = np.ascontiguousarray(inpt.PointData[self.inputArray]).astype('f4')

          ccode.SampleCartesian(self.samples, o, s, k, self.minSpacing, self.maxSpacing, 1, data, nRetained, rs, rd)

        else:
          inpt = vtk.vtkUnstructuredGrid.GetData(inInfoVec[0], 0)

          if inpt == None:
            print("Can only handle ImageData or UnstructuredGrid")
            return 1

          inpt = dsa.WrapDataObject(inpt)

          if np.min(inpt.CellTypes) < vtk.VTK_TETRA or np.max(inpt.CellTypes) >  vtk.VTK_TETRA:
            print("can handle only cartesian grids or unstructured grids containing only tetrahedra")
            return 1

          nCells  = inpt.GetNumberOfCells()
          points  = np.ascontiguousarray(inpt.Points).astype('f4')
          tets    = np.ascontiguousarray(inpt.Cells).astype('i4')
          data    = np.ascontiguousarray(inpt.PointData[self.inputArray]).astype('f4')

          ccode.SampleTetrahedra(self.samples, nCells, self.minSpacing, self.maxSpacing, points, tets, 1, data, nRetained, rs, rd)

        nSamples = ccode.GetNumberOfSamples()
        samples = np.zeros(nSamples * 3).astype('f4')
        ccode.GetSamples(samples)
        samples = samples.reshape(-1,3)

        print("Xreturn ", nSamples, " samples")

        so = dsa.WrapDataObject(vtk.vtkUnstructuredGrid())
        co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(nSamples).astype('i8')*2)
        ca = vtk.vtkCellArray()
        ca.SetCells(nSamples, dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([1]*nSamples, range(nSamples))).reshape((2*nSamples,)).astype('i8')))
        ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*nSamples).astype('u1'))
        so.VTKObject.SetCells(ct, co, ca)
        so.Points = dsa.numpy_support.numpy_to_vtk(samples, deep=1)

        print("hello")
        from vtk import vtkResampleWithDataSet
        r = vtkResampleWithDataSet()
        r.SetSourceData(inpt.VTKObject)
        r.SetInputData(so.VTKObject)
        r.Update()

        outpt = vtk.vtkUnstructuredGrid.GetData(outInfoVec, 0)
        outpt.ShallowCopy(r.GetOutput())

        if self.continuity:
          self.retained = r.GetOutput()

        del so
        del r

        return 1
        

