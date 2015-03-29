"""Analysis of ion binding.

"""

import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import mdsynthesis as mds

from MDAnalysis.core.distances import distance_array

# This function works, but it is very ugly, and has many repeating elements of
# code, making it difficult to maintain. To make it more maintainable we could
# define repeated elements as functions inside of it, but this would still be a
# bit sloppy since encapsulation is less clear due to scoping rules. 
# We can get more functionality out of a class.
def binding(universe, resids, ion, outfile='./iondata.h5', writeout=1000, start=0,
        stop=-1, step=1):
    """Get ion binding information for selected residues.

        :Arguments:
            *universe*
                universe to get data for
            *resids*
                resids to get ion binding information
            *ion*
                residue name of ion to get binding information for

        :Keywords:
            *outfile*
                filename for output data
            *writetout*
                how many frames to collect data for between writeouts
            *start*
                start frame
            *stop*
                end frame
            *step*
                frame step

    """
    
    ## pre-loop set up ###################
    # set up data structure
    data = dict()

    ions = universe.selectAtoms('resname {}'.format(ion))

    for resid in resids:
        data[resid] = dict()

        data[resid]['index'] = 0
        data[resid]['any_collected'] = False

        data[resid]['residue'] = universe.selectAtoms('resid {}'.format(resid))
    
        ## zero data structures
        data[resid]['times'] = np.zeros((writeout), dtype=float)

        # getting 4 closest ions for each frame
        data[resid]['data'] = np.zeros((writeout, 4), dtype=float)

        # atom index of first closest sodium, second, etc.
        data[resid]['ions'] = np.zeros((writeout, 4), dtype=int)
        data[resid]['any_collected'] = False

    ## main loop ###########################
    for ts in universe.trajectory[start:stop:step]:
        time = ts.time
        print "\rAnalyzing frame {}/{}, time {} ps".format(ts.frame, universe.trajectory.numframes, time),  
        for resid in data:

            # get data
            index = data[resid]['index']
            times = data[resid]['times']
            binding = data[resid]['data']
            ionarray = data[resid]['ions']

            residue = data[resid]['residue']

            times[index] = time
        
            dists = distance_array(ions.positions, residue.positions).min(axis=1)
            sort = dists.argsort()

            # get closest 4 ions
            binding[index, :] = dists[sort][:4]
            ionarray[index, :] = ions.indices()[sort][:4]

            # increment counter
            data[resid]['index'] = data[resid]['index'] + 1
            
            # make sure we write out at the end
            data[resid]['any_collected'] = True

            # if counter was equal to writeout, write it out!
            if data[resid]['index'] == writeout:

                ## writeout 
                times = data[resid]['times']
                binding = data[resid]['data']
                ionarray = data[resid]['ions']
                index = data[resid]['index']

                if index:
                    end = index
                # if index == 0 and any_collected == True, we've come full circle
                elif data[resid]['any_collected']:
                    end = -1
                # if any_collected is False, then we have an empty data structure; skip
                else:
                    return

                d = {
                     'dist1': binding[:end,0],
                     'ion1': ionarray[:end,0],
                     'dist2': binding[:end,1],
                     'ion2': ionarray[:end,1],
                     'dist3': binding[:end,2],
                     'ion3': ionarray[:end,2],
                     'dist4': binding[:end,3],
                     'ion4': ionarray[:end,3],
                    }

                cols = [
                        'dist1',
                        'ion1',
                        'dist2',
                        'ion2',
                        'dist3',
                        'ion3',
                        'dist4',
                        'ion4',
                       ]

                outdata = pd.DataFrame(d, columns=cols,
                            index=pd.Float64Index(times[:end], name='time'))

                # write to disk
                f = pd.HDFStore(outfile, mode='a')
                f.append(str(resid), outdata, data_columns=True)
                f.close()

                ## zero our existing arrays; saves re-allocation
                for item in ['times', 'data', 'ions']:
                    data[resid][item].fill(0)

                data[resid]['any_collected'] = False

                # reset counter
                data[resid]['index'] = 0
        

    ## post-loop wrap up ##################################
    for resid in data:
        ## writeout 
        times = data[resid]['times']
        binding = data[resid]['data']
        ionarray = data[resid]['ions']
        index = data[resid]['index']

        if index:
            end = index
        # if index == 0 and any_collected == True, we've come full circle
        elif data[resid]['any_collected']:
            end = -1
        # if any_collected is False, then we have an empty data structure; skip
        else:
            return

        d = {
             'dist1': binding[:end,0],
             'ion1': ionarray[:end,0],
             'dist2': binding[:end,1],
             'ion2': ionarray[:end,1],
             'dist3': binding[:end,2],
             'ion3': ionarray[:end,2],
             'dist4': binding[:end,3],
             'ion4': ionarray[:end,3],
            }

        cols = [
                'dist1',
                'ion1',
                'dist2',
                'ion2',
                'dist3',
                'ion3',
                'dist4',
                'ion4',
               ]

        outdata = pd.DataFrame(d, columns=cols,
                    index=pd.Float64Index(times[:end], name='time'))

        # write to disk
        f = pd.HDFStore(outfile, mode='a')
        f.append(str(resid), outdata, data_columns=True)
        f.close()

# The same anlaysis as above structured as a class. It happens to be slightly
# more code, but object attributes are easily distinguishable from method
# variables. It also allows us to use fancier elements such as properties (not
# used here)
class Binding(object):
    """Get ion binding information for selected residues.

    :WARNING: if there is more than one residue in your topology with the resid(s) used,
        then you can expect to see funny results

    Obtains the minimum distance (closest atoms) between the nearest 4 ions and the
    selected residues. The data is written as a pandas DataFrame to an HDF5
    file, with one DataFrame per residue. The DataFrames can be retrieved from
    the HDF5 file with the `pandas.HDFStore` object. Each DataFrame is indexed
    by simulation time (ps), and columns include distances of the nearest 4 ions
    and their atom ids for each frame.

    :NOTE: No sorting or duplicate removal is performed on the data stored. Each time
    the analysis is run, the results are simply appended to the existing data, if it exists.

    """
    def __init__(self, universe, resids, ion, outfile='./iondata.h5', writeout=1000):
        """Initialize analysis settings.

        :Arguments:
            *universe*
                universe to get data for
            *resids*
                resids to get ion binding information
            *ion*
                residue name of ion to get binding information for

        :Keywords:
            *outfile*
                filename for output data
            *writeout*
                how many frames to collect data for between writeouts

        """

        self.outfile = outfile

        self.universe = universe
        self.writeout = writeout
        self.ion = ion
        self.resids = resids

    def run(self, start=0, stop=-1, step=1):
        """Run analysis on universe.

        """
        ## pre-loop set up
        self.pre_loop()

        ## main loop
        for ts in self.universe.trajectory[start:stop:step]:
            self.time = ts.time
            print "\rAnalyzing frame {}/{}, time {} ps".format(ts.frame, self.universe.trajectory.numframes, self.time),  
            self.loop(ts)

        ## post-loop wrap up
        self.post_loop()

    def pre_loop(self):
        self.data = dict()
        self.ions = self.universe.selectAtoms('resname {}'.format(self.ion))

        for resid in self.resids:
            self.data[resid] = dict()

            self.data[resid]['index'] = 0
            self.data[resid]['any_collected'] = False

            self.data[resid]['residue'] = self.universe.selectAtoms('resid {}'.format(resid))
    
            # setup data structures
            self._setup(resid)

    def _setup(self, resid):
        """Set up data structures.

        It is expensive to allocate memory, so we create data structures for storing
        our values of fixed size and data type. We fill these up in the main loop,
        and when they are full, we package the data in a useful form (usually
        a pandas object), and write it out. We then re-zero our data structure.

        """
        self.data[resid]['times'] = np.zeros((self.writeout), dtype=float)

        # getting 4 closest ions for each frame
        self.data[resid]['data'] = np.zeros((self.writeout, 4), dtype=float)

        # atom index of first closest sodium, second, etc.
        self.data[resid]['ions'] = np.zeros((self.writeout, 4), dtype=int)
        self.data[resid]['any_collected'] = False

    def _zero(self, resid):
        """Zero data structures.

        It is expensive to allocate memory, so we create data structures for storing
        our values of fixed size and data type. We fill these up in the main loop,
        and when they are full, we package the data in a useful form (usually
        a pandas object), and write it out. We then re-zero our data structure.

        """
        # zero our existing arrays; saves re-allocation
        for item in ['times', 'data', 'ions']:
            self.data[resid][item].fill(0)

        self.data[resid]['any_collected'] = False

    def loop(self, ts, **kwargs):
        for resid in self.data:

            # get data
            self._ions(ts, resid)
            
            # make sure we write out at the end
            self.data[resid]['any_collected'] = True

            # if counter was equal to writeout, write it out!
            if self.data[resid]['index'] == self.writeout:
                self._writeout(resid)
                self._zero(resid)

                # reset counter
                self.data[resid]['index'] = 0

    def _writeout(self, resid):
        times = self.data[resid]['times']
        data = self.data[resid]['data']
        ions = self.data[resid]['ions']
        index = self.data[resid]['index']

        if index:
            end = index
        # if index == 0 and any_collected == True, we've come full circle
        elif self.data[resid]['any_collected']:
            end = -1
        # if any_collected is False, then we have an empty data structure; skip
        else:
            return

        d = {
             'dist1': data[:end,0],
             'ion1': ions[:end,0],
             'dist2': data[:end,1],
             'ion2': ions[:end,1],
             'dist3': data[:end,2],
             'ion3': ions[:end,2],
             'dist4': data[:end,3],
             'ion4': ions[:end,3],
            }

        cols = [
                'dist1',
                'ion1',
                'dist2',
                'ion2',
                'dist3',
                'ion3',
                'dist4',
                'ion4',
               ]

        outdata = pd.DataFrame(d, columns=cols,
                    index=pd.Float64Index(times[:end], name='time'))

        # write to disk
        f = pd.HDFStore(self.outfile, mode='a')
        f.append(str(resid), outdata, data_columns=True)
        f.close()
    
    def _ions(self, ts, resid):
        """Per resid calculations. Perhaps less efficient than doing all
           resids at once, but easier to write and more efficient when data
           already exists if checks are added for that. Plus, no scaling
           concerns.
    
        """
        index = self.data[resid]['index']
        times = self.data[resid]['times']
        data = self.data[resid]['data']
        ions = self.data[resid]['ions']

        residue = self.data[resid]['residue']

        times[index] = self.time
    
        dists = distance_array(self.ions.positions, residue.positions).min(axis=1)
        sort = dists.argsort()

        # get closest 4 ions
        data[index, :] = dists[sort][:4]
        ions[index, :] = self.ions.indices()[sort][:4]

        # increment counter
        self.data[resid]['index'] = self.data[resid]['index'] + 1

    def post_loop(self, **kwargs):
        for resid in self.data:
            self._writeout(resid)
