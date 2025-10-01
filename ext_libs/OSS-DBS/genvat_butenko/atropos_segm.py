import ants
import os
import sys

''' Atropos segmentation with possibility of multivariate input '''


def segment_with_atropos(imgPath, stimFolder):
    if type(imgPath) == list:
        img = []
        for img_path in imgPath:
            img.append(ants.image_read(img_path))
        mask = ants.get_mask(img[0])
    else:
        img = ants.image_read(imgPath)
        mask = ants.get_mask(img)
    segm_dict = ants.atropos(a=img, m='[0.2,1x1x1]', x=mask, v=1)
    segm_dict['segmentation'].image_write(os.path.join(stimFolder, 'segmask_raw_atropos.nii'))

    # we will re-number later in MATLAB


if __name__ == '__main__':
    # called from MATLAB
    # sys.argv[1] - path to the image
    # sys.argv[2] - stimulation folder

    if len(sys.argv) > 3:
        imgs = [sys.argv[1]]
        for img_i in range(3,len(sys.argv)):
            imgs.append(sys.argv[img_i])
            
        segment_with_atropos(imgs, sys.argv[2])
    else:
        segment_with_atropos(sys.argv[1], sys.argv[2])
