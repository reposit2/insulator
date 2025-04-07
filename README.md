## Demonstration: Prediction of Insulator-Associated DNA-Binding Proteins  

---  

This repository provides a demonstration of the prediction of insulator-associated DNA-binding proteins (DBPs), as described in our publication:

"Systematic discovery of directional regulatory motifs associated with human insulator sites."
[bioRxiv](https://doi.org/10.1101/2024.01.20.573595)

In this demonstration, we calculate DeepLIFT scores for DNA-binding sites (DBSs) of various DBPs in putative enhancers and promoters located upstream and downstream of each transcript. The method includes analysis of the well-known insulator-associated DBP, **CTCF**, as described in Figure 2d and the section *“Prediction of insulator-associated DBPs using the DeepLIFT tool”* in the Methods of our paper. The same approach is applied to other DBPs associated with insulator functions.

Due to file size limitations on GitHub, the required dataset is hosted on Zenodo.
[Zenodo repository](https://zenodo.org/record/8216164).  

---

### Setup  
  
#### System Environment  
- Ubuntu 20.04
- Python 3.8.12
- TensorFlow 2.6.3
- NVIDIA CUDA Toolkit 11.6
- GPU: GeForce RTX 3080 (10 GB memory)
- CPU: Intel Core i7-10700 (128 GB memory)

(*Note: The original paper used Python 3.8.10, TensorFlow 2.6.0, and CUDA Toolkit 11.4.2.*)

---

### Requirements

#### Data Files
Please download the following files from the [Zenodo repository](https://doi.org/10.5281/zenodo.8216164)

- `data.tgz`
- `mydata1.tgz`
- `mydata2.tgz`

Place all three files in the same directory as this repository, then extract them with the following command:

```bash
tar zxf *.tgz
```

#### Python Packages (tested versions)
- `numpy` (1.19.5)
- `pandas` (1.3.5)
- `scipy` (1.10.1)
- `deeplift` (0.6.12.0)
- `yaml` (6.0)
- `h5py` (3.1.0)
- `tensorflow` (2.6.3)

#### R Packages (tested versions)
- `tidyverse` (2.0.0)
- `ggpubr` (0.6.0)
- `broom` (1.0.4)

---

### Notes on DeepLIFT Compatibility
The original DeepLIFT package was developed for TensorFlow 1. For this project, we modified the codebase to support TensorFlow 2.
After installing DeepLIFT, we duplicated the original package directory, renamed it to deeplift_org, and applied our modifications within that directory.

```bash:
ls -l /home/user/.pyenv/versions/3.8.12/lib/python3.8/site-packages
drwxrwxr-x  6 user user   4096  3月 28 12:34 deeplift
drwxrwxr-x  2 user user   4096  3月 28 12:28 deeplift-0.6.13.0.dist-info
drwxrwxr-x  6 user user   4096  3月 28 12:34 deeplift_org
```

```diff
diff deeplift/conversion/kerasapi_conversion.py deeplift_org/conversion/kerasapi_conversion.py
386,387d385
+ <        if 'input' in layer_name:
+ <            continue
389,390c387,388
+ <            ("Layer "+str(layer_name)+" is in the layer names but not in the "
+ <             +" weights file which has layer names "+str(model_weights.keys()))
- >            ("Layer "+layer_name+" is in the layer names but not in the "
- >             +" weights file which has layer names "+model_weights.keys())
```

```diff
diff deeplift/layers/convolutional.py deeplift_org/layers/convolutional.py
8,9d7
+ < tf.compat.v1.disable_v2_behavior()
+ < tf.compat.v1.disable_eager_execution()
93c91
+ <         conv_without_bias = tf.compat.v1.nn.conv1d(
- >         conv_without_bias = tf.nn.conv1d(
253c251
+ <         conv_without_bias = tf.compat.v1.nn.conv2d(
- >         conv_without_bias = tf.nn.conv2d(
281c279
+ <                         tf.compat.v1.nn.conv2d_transpose(
- >                         tf.nn.conv2d_transpose(
288c286
+ <                        +tf.compat.v1.nn.conv2d_transpose(
- >                        +tf.nn.conv2d_transpose(
296c294
+ <                         tf.compat.v1.nn.conv2d_transpose(
- >                         tf.nn.conv2d_transpose(
303c301
+ <                        +tf.compat.v1.nn.conv2d_transpose(
- >                        +tf.nn.conv2d_transpose(
310c308
+ <             inp_mxts_increments += zero_inp_mask*tf.compat.v1.nn.conv2d_transpose(
- >             inp_mxts_increments += zero_inp_mask*tf.nn.conv2d_transpose(
```

```diff
diff deeplift/layers/core.py deeplift_org/layers/core.py
15,16c15
+ < tf.compat.v1.disable_v2_behavior()
+ < tf.compat.v1.disable_eager_execution() 
66c65
+ <         self._pos_mxts = tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >         self._pos_mxts = tf.zeros_like(tensor=self.get_activation_vars(),
68c67
+ <         self._neg_mxts = tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >         self._neg_mxts = tf.zeros_like(tensor=self.get_activation_vars(),
208c207
+ <         self._activation_vars = tf.compat.v1.placeholder(
- >         self._activation_vars = tf.placeholder(
216c215
+ <         return tf.compat.v1.placeholder(dtype=tf.float32,
- >         return tf.placeholder(dtype=tf.float32,
456c455
+ <              tf.compat.v1.variables_initializer([self.task_vector])) 
- >              tf.variables_initializer([self.task_vector])) 
472c471
+ <         task_vector_update = tf.compat.v1.assign(self.task_vector,
- >         task_vector_update = tf.assign(self.task_vector,
474c473
+ <         task_vector_update = tf.compat.v1.scatter_update(
- >         task_vector_update = tf.scatter_update(
```

```diff
diff deeplift/layers/helper_functions.py deeplift_org/layers/helper_functions.py
2,4d1
+ < tf.compat.v1.disable_v2_behavior()
+ < tf.compat.v1.disable_eager_execution()
24c21
+ <     return tf.squeeze(tf.compat.v1.nn.conv2d_transpose(
- >     return tf.squeeze(tf.nn.conv2d_transpose(
```

```diff
diff deeplift/layers/pooling.py deeplift_org/layers/pooling.py
10,11c10
+ < tf.compat.v1.disable_v2_behavior()
+ < tf.compat.v1.disable_eager_execution()
64c63
+ <                 tf.compat.v1.nn.max_pool(value=tf.expand_dims(input_act_vars,1),
- >                 tf.nn.max_pool(value=tf.expand_dims(input_act_vars,1),
74c73
+ <         return tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >         return tf.zeros_like(tensor=self.get_activation_vars(),
76c75
+ <                tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >                tf.zeros_like(tensor=self.get_activation_vars(),
117c116
+ <         return tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >         return tf.zeros_like(tensor=self.get_activation_vars(),
119c118
+ <                tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >                tf.zeros_like(tensor=self.get_activation_vars(),
196,197c195
+ < #        width = self._get_input_activation_vars().get_shape().as_list()[1]
+ <         width = self._get_input_activation_vars().get_shape().as_list()[2]
- >         width = self._get_input_activation_vars().get_shape().as_list()[1]
266c264
+ <         to_return = tf.compat.v1.nn.max_pool(value=input_act_vars,
- >         to_return = tf.nn.max_pool(value=input_act_vars,
280c278
+ <         return tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >         return tf.zeros_like(tensor=self.get_activation_vars(),
282c280
+ <                tf.compat.v1.zeros_like(tensor=self.get_activation_vars(),
- >                tf.zeros_like(tensor=self.get_activation_vars(),
```

---

### Run
  
```
bash calc_deepliftscore.sh
bash stat_deepliftscore.sh
perl stat2_deepliftscore.pl 1
perl stat2_deepliftscore.pl 2
```

---

### Outputs

We analyzed two sets of DNA-binding sites (DBSs) for CTCF, based on the following selection criteria:

- **DBS Selection Criterion S1**
- **DBS Selection Criterion S2**

These criteria are defined and illustrated in **Supplementary Fig. 1** of our paper.
The **absolute z-scores** for DBSs selected using **Criterion S2** are shown in **Figure 3d**.

```
perl stat2_deepliftscore.pl 1

CTCF	FR_vs._RF	FR_vs._FF	FR_vs.RR	FR_vs.X
zscore	30.5086444636761	22.6982073925036	22.1370605812343	-4.15304647173935
pvalue	 2.001291193103142e-204	 4.666373842989266e-114	 1.390119461104552e-108	 3.280780828751686e-05
mean_difference	0.000129712500594734	0.000106661611432207	0.000122422833153355	1.53591758089962e-05
mean	0.0024391885341362956	0.0024119834276677098	0.0024181798004967595	0.0023663710668580524	0.002309476033541562	0.002305321816235503	0.002295756967343405	0.0023510118910490562
rvalue	0.02175830785424293	0.02137963064522097	0.02080993035315193	0.004149985121522093
#_scores	1352987	766181	769080 	648057	1270983	737247	740290	688467

perl stat2_deepliftscore.pl 2

CTCF	FR_vs._RF	FR_vs._FF	FR_vs.RR	FR_vs.X
zscore	29.3514936805777	20.0532840167669	20.0667097837562	-10.908618128604
pvalue	 2.286905905408163e-189	 1.88951021507238e-89	 1.442432336817605e-89	 1.0483782537220221e-27
mean_difference	0.0001211693880644	9.09719290027194e-05	0.000110929932730617	4.12370877662651e-05
mean	0.0023831139855248745	0.0023536584654016545	0.0023575815683481555	0.0023206333968334067	0.002261944597460475	0.002262686536398935	0.0022466516356175383	0.0022793963090671416
rvalue	0.02214496692055956	0.02006341963250835	0.020128639347794408	0.011655044924477997
#_scores	1211945	684926 	679036 	513268	1133046	648080	646886	677371
```

---

### Reference and License

#### DEcode
Tasaki, S., Gaiteri, C., Mostafavi, S., and Wang, Y.  
**Deep learning decodes the principles of differential gene expression.**  
*Nature Machine Intelligence* (2020).  
[Link to paper](https://doi.org/10.1038/s42256-020-0201-6)  
License: BSD 3-Clause

#### DeepLIFT
Shrikumar, A., Greenside, P., and Kundaje, A.  
**Learning Important Features Through Propagating Activation Differences.**  
*arXiv* (2017).  
[Link to paper](https://arxiv.org/abs/1704.02685)
License: MIT



