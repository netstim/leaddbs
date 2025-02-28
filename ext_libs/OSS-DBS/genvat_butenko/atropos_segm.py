import ants
import os
import sys

''' Atropos segmentation, for now with one image '''


def segment_with_atropos(imgPath, stimFolder):
    img = ants.image_read(imgPath)
    mask = ants.get_mask(img)
    segm_dict = ants.atropos(a=img, m='[0.2,1x1x1]', x=mask, v=1)
    segm_dict['segmentation'].image_write(os.path.join(stimFolder, 'segmask_raw_atropos.nii'))

    # we will re-number later in MATLAB


if __name__ == '__main__':
    # called from MATLAB
    # sys.argv[1] - path to the image
    # sys.argv[2] - stimulation folder

    segment_with_atropos(sys.argv[1], sys.argv[2])
