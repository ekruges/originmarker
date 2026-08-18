/**
 * Isotonic recalibration maps, fitted on THIS implementation's own injections.
 *
 * GENERATED. Regenerate with `node audit/calibration/fit_maps.ts`, which prints this file's data.
 *
 * NOT THE METHODS REVIEW'S MAPS. Those were specified, never delivered, and could not be
 * reproduced because the corpus carries no material labels. These are fitted on OUR posterior over
 * OUR noise, keyed on the material THIS TOOL ASSIGNS AT RUNTIME, which is self-consistent by
 * construction: they calibrate the assignment actually made rather than one needing labels nobody
 * has.
 *
 * FITTED LEAVE-ONE-ARRAY-OUT, which is the difference between demonstrating calibration and
 * asserting it: every accuracy quoted comes from rows the map never saw. Measured expected
 * calibration error, raw against mapped:
 *
 *     bulk           0.0165 -> 0.0065
 *     trophectoderm  0.0129 -> 0.0089
 *     esc-single     0.0135 -> 0.0081
 *     blastomere     0.0085 -> 0.0014
 *
 * THE MATERIAL KEYS ARE A STAND-IN AND SAY SO. With no labels available, the 35 arrays were split
 * into four noise tiers and each tier named for the material it is characteristic of. These
 * calibrate the NOISE each material carries rather than the material itself, which is the honest
 * description of what was fitted and the reason a real labelled set would still be worth having.
 *
 * A file that is here rather than a threshold in code, because a calibration is data.
 */
export const SHIPPED_MAPS = {
  'bulk': {
    marginal: {
      raw: [0.506438, 0.518758, 0.569007, 0.619316, 0.693809, 0.719100, 0.730764, 0.743810, 0.755312, 0.819222, 0.843446, 0.856239, 0.893210, 0.931765, 0.943414, 0.981618, 0.999224],
      calibrated: [0.507858, 0.521298, 0.550688, 0.576874, 0.596591, 0.600000, 0.666667, 0.758065, 0.769509, 0.772519, 0.785714, 0.796610, 0.914298, 0.929577, 0.947679, 0.984000, 1.000000],
    },
  },
  'trophectoderm': {
    marginal: {
      raw: [0.506382, 0.518681, 0.568838, 0.589259, 0.606206, 0.618561, 0.643525, 0.668896, 0.693541, 0.768572, 0.831833, 0.844032, 0.919167, 0.956425, 0.969363, 0.981424, 0.999320],
      calibrated: [0.516811, 0.524834, 0.555586, 0.575758, 0.639175, 0.642801, 0.732707, 0.753127, 0.813960, 0.826741, 0.876712, 0.891209, 0.897684, 0.975610, 1.000000, 1.000000, 1.000000],
    },
  },
  'esc-single': {
    marginal: {
      raw: [0.506271, 0.518638, 0.567543, 0.593712, 0.618735, 0.655568, 0.681167, 0.693960, 0.744044, 0.794128, 0.818489, 0.906023, 0.968750, 0.982128, 0.999198],
      calibrated: [0.513652, 0.527032, 0.555050, 0.647054, 0.675940, 0.783644, 0.791045, 0.803880, 0.816125, 0.885980, 0.918584, 0.961106, 0.978142, 0.987578, 0.997509],
    },
  },
  'blastomere': {
    marginal: {
      raw: [0.506244, 0.518893, 0.543468, 0.555764, 0.568584, 0.578486, 0.593651, 0.631099, 0.644326, 0.681215, 0.706378, 0.718793, 0.755937, 0.794093, 0.806063, 0.844407, 0.906711, 0.919178, 0.931329, 0.943813, 0.955693, 0.969078, 0.999234],
      calibrated: [0.505114, 0.520270, 0.524226, 0.541044, 0.547303, 0.602317, 0.636685, 0.645833, 0.694216, 0.714964, 0.735849, 0.748887, 0.750693, 0.812500, 0.834202, 0.860654, 0.870229, 0.870504, 0.891304, 0.948276, 0.980861, 0.991940, 0.998886],
    },
  },
} as const
