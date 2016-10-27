
###############################

How to Build the Model:

Run survivalFactorizationEM.SurvivalFactorizationEM_Runner providing a path of a configuration file. Given a dataset, this script generates several instances of class survivalFactorizationEM.SurvivalFactorizationEM_Model.
The configuration file must to be written according to the .”properties” syntax. The fields are:

	n_factors = <list of the number of latent factors (topics)>
	output = <output path prefix>
	max_iterations = <maximum number of iterations>
	assignment_file = <topic file>
	event_file = <cascades>
	[ content_file = <texts> ]

where

<list of the number of latent factors (topics)>: an integer list separated by “;”
	e.g. n_factors = 2;4;8;16;32;64;128
This list sets the number of models to build, one for each number of factors.

<output path prefix>: a String
	e.g. output = resources/datasets/synth/models/Synth
This string contains two elements:
	- the path of folder where the built models (in the example “resources/datasets/synth/models”) will be stored
	- the prefix of the name of the file which will contain the model (in the example “Synth”).
For each topic number in <list of the number of latent factors (topics)>, a model, with the corresponding number of factors, will be created and will be stored in the folder; the name of the model file is a concatenation of the prefix + “_” + <number of factors> + “.model”

<maximum number of iterations>: an integer
	e.g. max_iterations = 1000
The fix point iterations will continue until convergence or when this number (burn-in phase included) is reached

<topic file>: a String
	e.g. assignment_file = resources/datasets/synth/models/Synth
This field is similar to <output path prefix>. Each assignment file will contain the association cascade - topic for each cascade

<cascades>: a String
	e.g. event_file = resources/datasets/Synth/cascades_training.txt
The name of a file containing the cascades of events (e.g. tweets) exploited to build the model

<texts>: a String
	e.g. event_file = resources/datasets/Synth/text_training.txt
The name of a file containing the text information for each cascade. Note: this field is optional


###############################

Cascade file format. A text document containing this information:

NodeId	CascadeId	TimeStamp
1449	1	1222254982000
6930	1	1222277866000
2466	2	1222281238000
…

Each row is an activation, the separator is the tab character “\t”


###############################

Text file format. A text document containing this information:

WordId	CascadeId	Frequency
1248	1	1
804	1	5
6788	2	3
8134	2	1
…

Each row is an word in an event, the separator is the tab character “\t”


###############################

How to perform the Network Reconstruction:

Run survivalFactorizationEM.FullTestNetworkReconstruction providing a path of a configuration file. Given a folder containing built models and a test network, this script will generate the network reconstruction for each model. The configuration file must to be written according to the .”properties” syntax. The fields are:

model_folder = <model folder>
model_files = <models>
test_file = <test network>
output_folder = <output folder>
output_files = <output files>

where:

<model folder>: a String
	e.g. model_folder = resources/datasets/synth/models
The folder containing the built models

<models>: a String list whose elements are separated by “;”
	e.g. model_files = Synth_2f.model;Synth_4f.model;Synth_8f.model;Synth_16f.model;Synth_32f.model;Synth_64f.model;Synth_128f.model
This list contains the file names where the models are stored

<test network>: a String
	e.g. test_file = resources/datasets/Synth/s2/links_reduced-FF1400.remapped_two_hops
The file containing the test network to reconstruct. Note: the syntax of the test file is equal to the Cascade file format.

<output folder>: a String
	e.g. output_folder = resources/datasets/Synth/preds
The folder where to put the reconstructed network files (one for each model)

<output files>: a String list whose elements are separated by “;”
	e.g. output_files = Synth_2f.pred;Synth_4f.pred;Synth_8f.pred;Synth_16f.pred;Synth_32f.pred;Synth_64f.pred;Synth_128f.pred
The file names of the reconstructed networks


###############################

Network reconstruction output prediction file format.

Prediction	ActualClass
1.15421803E-6	1
2.27729428E-6	1
1.16779013E-7	2
…

Each row is a prediction, the separator is the tab character “\t”