# single_patch_pair.py
# Utility functions for generating lesion-based patch pairs from MRI slices

import numpy as np
import scipy

# Internal mask cache
MASK_CACHE = {}

def load_kspace_masks(mask_dir="./Mask"):
    """
    Loads and caches 2x and 4x undersampling masks from .npy files.
    Must be called once before using apply_mask_k().
    """
    global MASK_CACHE
    if "2x" not in MASK_CACHE:
        m2x = np.load(f"{mask_dir}/mask_2x_cartesian_for_260x311.npy")
        m2x = scipy.ndimage.shift(m2x, [m2x.shape[0]//2, m2x.shape[1]//2], mode='wrap')
        MASK_CACHE["2x"] = m2x
    if "4x" not in MASK_CACHE:
        m4x = np.load(f"{mask_dir}/mask_4x_cartesian_for_260x311.npy")
        m4x = scipy.ndimage.shift(m4x, [m4x.shape[0]//2, m4x.shape[1]//2], mode='wrap')
        MASK_CACHE["4x"] = m4x


def mri_valid_regions(image, mode='both'):
    lower = 3000.
    upper = 5000.
    if (mode == 'left'):
        upper = 4000.
    elif (mode == 'right'):
        lower = 4000.
    mod_image = np.where((image >= lower) & (image < upper), 1, 0).astype(np.int8)
    return mod_image


def find_spots(img, kernel, num=2, unpad=20):
    """
    Returns a list of [y, x] points where a lesion represented
    by a 64x64 kernel may be placed in the image. This guarantees
    all nonzero parts of the lesion go over ones in the image.

    Parameters:
        img (ndarray):The image where patches are to be found. This only
            works if the image is composed of ones and zeros.
        kernel (ndarray):Kernel containing a lesion to be placed on the
            image. Must be 64x64 specifically.
        num (int):Number of points to find. Default is 2.
        unpad (int):Number of pixels on each side to exclude from the
            search for faster execution. Default is 20.

    Returns:
        retpts (list):List of 2-entry lists, each containing the [y,x]
            point of the bottom right pixel where the kernel may be overlaid
            on the image. The slice may be removed as [y-64+1:y+1,x-64+1:x+1].
    """
    kcopy = np.where(kernel > 0.1 * kernel.max(), 1, 0)
    retpts = []
    if not (unpad):
        unpad = 0
    for a in range(num):
        flag = 0
        spincount = 0
        while (flag == 0):
            prop_y = np.random.randint(64 + unpad, img.shape[0] + 1 - unpad)
            prop_x = np.random.randint(64 + unpad, img.shape[1] + 1 - unpad)
            # if ((prop_y > 63) and (prop_y <= img.shape[0] + 1) and
            #    (prop_x > 63) and (prop_x <= img.shape[1] + 1)):
            if (np.multiply(img[prop_y - 63:prop_y + 1, prop_x - 63:prop_x + 1], kcopy).sum() == kcopy.sum()):
                retpts.append([prop_y, prop_x])
                flag = 1
            spincount = spincount + 1
            if (spincount >= 5000):
                raise ValueError(
                    "Couldn't find a point in 1000 attempts: make sure the image is binary and has patches large enough for the lesion")
    if (flag == 0):
        print("Error! Didn't find enough points!")
    return retpts


def add_lesion(img, kernel, loc, copy=False):
    """
    The values in the kernel are added to the image in-place, unless copy
    is set to True.

    add_lesion(img, kernel, loc, copy=False)
    Parameters:
        img (ndarray):The image where a lesion will be added.
        kernel (ndarray):Kernel containing a lesion to be placed on the
            image. Must be 64x64 specifically.
        loc (list of two ints):Location of one less in each axis than
            the lower-right corner of where the kernel will be added.
        copy (boolean): Whether the kernel will be added to a copy of the
            image rather than in-place. Default is False.

    Returns:
        Resulting image with lesion/signal added. No need to use this
        if copy is False.
    """
    if (copy == True):
        ret_img = img.copy()
    else:
        ret_img = img.view()
    ret_img[loc[0] - 64 + 1:loc[0] + 1, loc[1] - 64 + 1:loc[1] + 1] = ret_img[loc[0] - 64 + 1:loc[0] + 1,
                                                                      loc[1] - 64 + 1:loc[1] + 1] + kernel
    return ret_img


# Load sampling masks
kspace__mask2x = np.load("Mask/mask_2x_cartesian_for_260x311.npy")
kspace__mask2x = scipy.ndimage.shift(kspace__mask2x, [kspace__mask2x.shape[0] // 2, kspace__mask2x.shape[1] // 2],
                                     mode='wrap')
kspace__mask4x = np.load("Mask/mask_4x_cartesian_for_260x311.npy")
kspace__mask4x = scipy.ndimage.shift(kspace__mask4x, [kspace__mask4x.shape[0] // 2, kspace__mask4x.shape[1] // 2],
                                     mode='wrap')


def forward_k(img):
    """
    Wrapper for scipy.fft.fft2, cannot be done in-place.
    """
    return scipy.fft.fft2(img)


def add_noise_k(kspace, std):
    """
    Adds zero-mean i.i.d. Gaussian noise with a specified standard
    deviation.
    """
    fxn_noise = np.random.normal(0, std, kspace.shape)
    return (kspace + fxn_noise)


def apply_mask_k(kspace, mask=2):
    """
    Multiplies by sampling mask.
    """
    if ((mask == 2) or (mask == '2')):
        return np.multiply(kspace, kspace__mask2x)
    elif ((mask == 4) or (mask == '4')):
        return np.multiply(kspace, kspace__mask4x)
    else:
        raise ValueError("Invalid mask: set mask=2 or mask=4")


def inverse_k(kspace):
    """
    Wrapper for scipy.fft.ifft2, cannot be done in-place.
    """
    return scipy.fft.ifft2(kspace)


indices_shifted = np.indices((128, 128)) - 63.5


def get_lesion_kernel(rot=30, intensity=1):
    """
    get_lesion_kernel(rot=30)
        Returns a 64x64 kernel with a elliptical lesion with standard
        deviation of 4 along x and 2 along y, rotated by angle rot (degrees).
        The peak magnitude of the Gaussian ellipses is 1 (not mulitiplied by
        a coefficient). Default rotation angle is 30 degrees.
        Gaussian ellipses is calculated based on eq. 18 from
            https://pmc.ncbi.nlm.nih.gov/articles/PMC10520791/#sec4
    """
    indices_sq_rotx = np.square(scipy.ndimage.rotate(indices_shifted[1, :, :], rot, reshape=False))
    indices_sq_rotx = indices_sq_rotx[32:96, 32:96]
    indices_sq_roty = np.square(scipy.ndimage.rotate(indices_shifted[0, :, :], rot, reshape=False))
    indices_sq_roty = indices_sq_roty[32:96, 32:96]

    sigma_x, sigma_y = 4, 2
    gauss = np.exp((indices_sq_rotx[:, :] / (-2 * (sigma_x ** 2))) +
                   (indices_sq_roty[:, :] / (-2 * (sigma_y ** 2))))
    return gauss * intensity

def generate_single_patch_pair(MRI_imgs, MRI_seg, i=0, kernel=None, noise=1e-3, aceleration=2):
    """
    Generate a signal-present and signal-absent patch pair from one MRI image.
    Returns:
        signal_patch (64x64), nonsignal_patch (64x64), metadata (i, x, y), reference_patch (64x64)
    """
    img = MRI_imgs[i, :, :].copy()
    seg = MRI_seg[i, :, :].copy()

    # 1. Get binary white matter mask
    valid_mask = mri_valid_regions(seg)

    # 2. Get lesion placement locations
    lesion_spots = find_spots(valid_mask, kernel, 2) # Find to lesion locations
    if len(lesion_spots) < 2:
        print("Not enough valid lesion spots.")
        return None

    lesion_loc = lesion_spots[0] # First possible location
    alt_loc = lesion_spots[1] # Second possible location

    # 3. Preserve clean copy of image ( for non-signal patch)
    original_img = img.copy()

    # 4. signal-present patch, forward_k and inverse_k already have shift
    add_lesion(img, kernel, lesion_loc)
    img_k = forward_k(img)
    img_k.real = add_noise_k(img_k.real, noise * img_k.real.max())
    img_k.imag = add_noise_k(img_k.imag, noise * img_k.imag.max())
    img_k = apply_mask_k(img_k, aceleration)
    recon_array = inverse_k(img_k).real
    y, x = lesion_loc
    signal_patch = recon_array[y-63:y+1, x-63:x+1]

    # === signal-absent patch ===
    orig_k = forward_k(original_img)
    orig_k.real = add_noise_k(orig_k.real, noise * orig_k.real.max())
    orig_k.imag = add_noise_k(orig_k.imag, noise * orig_k.imag.max())
    orig_k = apply_mask_k(orig_k, aceleration)
    orig_recon = inverse_k(orig_k).real
    y_ns, x_ns = alt_loc
    nonsignal_patch = orig_recon[y_ns-63:y_ns+1, x_ns-63:x_ns+1]

    # reference patch from original image (no lesion, no noise, no downsampling)
    reference_patch = original_img[y_ns-63:y_ns+1, x_ns-63:x_ns+1]

    return signal_patch, nonsignal_patch, (i, x, y), reference_patch