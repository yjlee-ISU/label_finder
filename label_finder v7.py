# in vivo isotope labling finder. Lee Lab Jan 2022
#
# Read Exported file from QualBrowser: copy exact mass list from QualBrowser and save as a tab delimited text file
#      It has four columns: m/z, intensity, relative intensity, and noise.
# Revised by Young Jin Lee on Jan 12, 2022
# Further revised by Young Jin Lee on Dec 1, 2022 to add comments and explanation.
# Debugged remaining issues. Changed to accept file names from conol. Dec 13, 2022.
# Changed the filtering tolerance to 2ppm for all signals, and 1.5 ppm for S/N > 30. Dec 19, 2022
#
# It is written for 13C-labeling but can be modified as necessary for other labeling.
#
# Overall procedure:
# 1. Read both 12C & 13C data sets. They are ascii data exported from QualBrowser. Modify if there is formatting difference.
# 2. Remove 13C1 label in 12C data. This is to avoid duplicate labeling; i.e., 13C2 peak in 13C data is also assigned as M+1 of 13C1 peak.
# 3. Remove common peaks in 13C data that are also present in 12C data. This is to remove contamination peaks that are not biological but it will also remove monoisotope peaks.
# 4. Assign each 13C labeled peak in 13C data to the peaks in 12C data for up to the max # of labeling if within tolerance (default = 3ppm).
# 5. Add back monoisotope peak in 13C data if they are removed in process 3.
# 6. Remove the labeling if there is only one label or if there are two labels but monoisotope and 13C1 label.
# i.e., We accept only if there is at least two labeling exist and the max labeling is not one.

def input_data(file1):
    spectra, i, factor, blank_line = [], 0, 1, 0

    f=open(file1,'r')

# Remove the first few dummy lines in the QualBrowser exported dataset that are none-data.
    i = 0
    print('Number of blank (or dummy) lines: ')
    blank_line = int(input())
    while i < blank_line:
        f.readline()
        i += 1

# Read each line and store 'm/z' and 'intensity' as array only if they are greather than the noise * factor. 
# Factor is 1 by default. It can be increased if you want to increase signal to noise cutoff. 
    for line in f:
        mz, intensity, relative, noise = line.split()
        mz, intensity, relative, noise  = float(mz), float(intensity), float(relative), float(noise)
        if noise < 0.5:
                noise = 0.5
        if intensity > noise:
            spectra.append((mz,intensity, noise))

    spectra.sort()
    f.close()

    return spectra

# Main

import sys
import math

# 'tolerance' is mass tolerance used for deleting masses from labeled data that appears in control, and for assigning labels e.g. 0.0000025 is 2.5ppm.
tolerance, spectra_peak, spectra_bkg = 0.000002, [], []

# Step1: Read non-labeled file, 12C data ('file0', tab delimited text file), and store as two-dimensional array 'spectra_bkg'.
print ('Filename for C12 data: ')
file0 = input()
spectra_bkg = input_data(file0)
num_bkg = len(spectra_bkg)

# Read labeled file, 13C data ('file1', tab delimited text file), and store as two-dimensional array 'spectra_peak'.
# 'spectra_peak0' is a back-up array needed to restore monoisotope peaks that are removed.
print ('Filename for C13 data: ')
file1 = input()
spectra_peak = input_data(file1)

num_peak = len(spectra_peak)

print ('Output Filename: ')
file_out = input()
f2 = open(file_out,'w')

# Change 'label_mass' depending on isotope labeling. e.g. 13C = 1.0033548, D = 1.006277, 15N = 0.997035

# Step 2: Remove labeled peaks within unlabeled data (e.g., 13C natural isotope peaks in 12C data).
# Currently it removes only 13C1 (count < 2) as it seems to be sufficient but can be extended for 13C2 or higher if necessary.

j, k, count, label_mass = 0, 1, 1, 1.0033548

while ((j+k) < num_bkg):
    mz0 = spectra_bkg[j][0]
    while ((j+k) < num_bkg) and (count < 2):
        mz1 = spectra_bkg[j+k][0]
        mz_label = mz0 + count * label_mass

        mz_low = mz_label*(1-tolerance)
        mz_high = mz_label*(1+tolerance)
        if (mz1 < mz_low):
            k += 1
        elif (mz1 < mz_high):
            spectra_bkg.pop(j+k)
            num_bkg -= 1
        else:
            count +=1
    j += 1
    k = 1
    count = 1

# Step 3: Remove common peaks between labeled and unlabeled in labeled data set. Note: mono-isotope peaks will be also lost.
#         If labeled peak is larger than twice of unlabeled, keep the peaks in the assumption it may not be a contamination.

num_bkg = len(spectra_bkg)

i, j, removed = 0, 0, []
while (i < num_peak) and (j < num_bkg):
    mz_peak = spectra_peak[i][0]
    mass_tol = mz_peak*tolerance
    mass_low = mz_peak - mass_tol
    mass_high = mz_peak + mass_tol
    mz_bkg = spectra_bkg[j][0]

    peak_int = spectra_peak[i][1]
    bkg_int = spectra_bkg[j][1]
    if mz_bkg > mass_low:
        if mz_bkg > mass_high:
            i += 1
        else:
            if peak_int < (2*bkg_int):
                removed.append(spectra_peak.pop(i))
                num_peak -= 1
            i += 1
    else:
        j += 1

# Step 4: Main part for assigning labeled to unlabeled.
# for each j in unlabeled, find all i in labeled with count number of labels
# Adjust 'effective_atom' as necessary. This is to determine an upperlimit of possible labeling: e.g., 14 for 13C labeling (one carbon per every CH2).
#   7 (=14/2) for 2H labeling (two hydrogens per every CH2).

i, j, i0, effective_atom = 0, 0, 0, 14.0
final = []
num_peak = len(spectra_peak)

while (j < num_bkg):
    count = 0
    mz_bkg = spectra_bkg[j][0]
    int_bkg = spectra_bkg[j][1]
    mz_tol1 = mz_bkg * tolerance
    mz_tol2 = mz_bkg * tolerance * 0.75
    while ((i < num_peak) and (count < (mz_bkg/effective_atom))):
        mz_peak0 = spectra_peak[i][0]
        mz_int = spectra_peak[i][1]
        mz_peak = mz_bkg + label_mass * count

        bkg_low = mz_peak - mz_tol1
        bkg_high = mz_peak + mz_tol1
        bkg_low2 = mz_peak - mz_tol2
        bkg_high2 = mz_peak + mz_tol2

        if (mz_peak0 < bkg_low):
            if count == 0:
                i0 = i
            i += 1
        elif (mz_peak0 < bkg_high):
            noise = spectra_peak[i][2]
            if mz_int > 30*noise:
                if (mz_peak0 > bkg_low2) and (mz_peak0 < bkg_high2):
                    final.append((mz_bkg, mz_peak0, count, mz_int, int_bkg))               
            else:
                final.append((mz_bkg, mz_peak0, count, mz_int, int_bkg))
            i += 1
            count += 1
        else:
            count += 1
    j += 1
    i = i0
    count = 0

# Step 5: Put back monoisotope peaks to 13C data if present in back-up labeled data.

i, j = 0, 0
num_assigned = len(final)
num_removed = len(removed)

while (j < num_assigned) and (i < num_removed):
    unlabeled = final[j][0]
    unlabel_low = unlabeled*(1-tolerance)
    unlabel_high = unlabeled*(1+tolerance)
    labeled = removed[i][0]
    count = final[j][2]

    if labeled > unlabel_low:
        if labeled > unlabel_high:
            j += 1
        else:
            if count != 0:
                labeled_intensity = removed[i][1]
                unlabeled_intensity = final[j][4]
                final.append((unlabeled, labeled, 0, labeled_intensity, unlabeled_intensity))
            i += 1
            j += 1
    else:
        i += 1

final.sort()

# Step 6: Remove the labeling if there is only one label or if there is two labels but 13C0 and 13C1 label.

i, k = 0, 1
num_assigned = len(final)

while (i+k) < num_assigned:
    current = final[i][0]
    next_peak = final[i+k][0]
    if (current != next_peak):
        if k == 1:
            final.pop(i)
            num_assigned -= 1
        elif (k == 2) and (final[i+1][2] == 1):
            final.pop(i)
            final.pop(i)
            num_assigned -= 2
            k = 1
        else:
            i += k
            k = 1
    else:
        k += 1
if final[i-1][0] != final[i][0]:
    final.pop(i)

# Save files

set(final)
i = 0
num_assigned = len(final)
s = 'Unlabeled' + '\t' + 'Labeled' +'\t'+'# of Labels'+ '\t' + 'Labeled intensity' + '\t' + 'Unlabeled intensity' +'\n'
f2.writelines(s)
while i < num_assigned:
    s =str(final[i][0])  + '\t' + str(final[i][1]) + '\t'+str(final[i][2])+ '\t' + str(final[i][3])+'\t' + str(final[i][4])+'\n'
    f2.writelines(s)
    i += 1
f2.close()
