#ifndef PSY_CSF_WEIGHTS_H
#define PSY_CSF_WEIGHTS_H

/* Initialize CSF weight tables. Called automatically on first use. */
void psy_csf_init(void);

/* Get CSF weight for a specific coefficient position in a block.
 * size: 4 or 8. coeff_idx: row*size+col in the coefficient block. */
float psy_csf_weight(int size, int coeff_idx);

/* Get pointer to the full CSF weight table for a block size. */
const float *psy_csf_table(int size);

#endif
