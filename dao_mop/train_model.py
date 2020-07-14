"""
Train up a CNN using the Keras API to TensorFlow. Also makes use of sklearn to build the test / validation set.

Keras is TensorFlow's API for building and training deep learning models.
It's used for fast prototyping, state-of-the-art research, and production.

TensorFlow is an open source package for the machine laerning.
TF has a lot of advantages compared to other ML packages.
For example, TF automatically calculates the derivates, and it express what it is doing well, and it can use CPU, GPU
and TPU. We can use it with Python, C++, Jana or R.

Scikit learn is a Python module for machine learning built on SciPy.
1. Supervised Learning module, and Unsupervised Learning module
2. model selection and validation module
3. Data conversion and Data loading module
"""
import argparse
import logging
import sys
import numpy as np
import keras
from keras.layers import BatchNormalization
from keras.layers import Dense, Dropout
from keras.layers import Flatten
from keras.layers import Input
from keras.layers.convolutional import Conv2D
from keras.layers.pooling import MaxPooling2D
from keras.models import Model
from matplotlib import pyplot as plt
from sklearn import model_selection as md
from . import data_model

CUTOUT_DIMENSION = 64
NUM_CHANNELS = 2


def plot_training_outcome(history, output_file_base=None):
    """
    Plot the accurcy and loss values as a function of epoch for the given model training history.
    :param history: Output of Model.fit
    :param output_file_base: name of file to save plots to (_accuracy.pdf and _loss.pdf is appended) if None then
    display plot with 'show' rather than save to file.
    :return: None
    """

    # Plot the 'accuracy' (true positive) rate for both training and validation.
    plt.plot(history['accuracy'])
    plt.plot(history['val_accuracy'])
    plt.title('model accuracy with Aram\'s data')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    if output_file_base is None:
        plt.show()
    else:
        plt.savefig("{}_accuracy.pdf".format(output_file_base))

    # Plot the 'loss' (false negative) rate for both training and validation.
    plt.plot(history['loss'])
    plt.plot(history['val_loss'])
    plt.title('model loss with Aram\'s data')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    if output_file_base is None:
        plt.show()
    else:
        plt.savefig("{}_loss.pdf".format(output_file_base))


class TrainingPlot(keras.callbacks.Callback):

    def __init__(self):
        self.losses = []
        self.acc = []
        self.val_losses = []
        self.val_acc = []
        self.logs = []
        self.x = []
        self.i = 1
        super().__init__()

    # This function is called when the training begins
    def on_train_begin(self, logs=None):
        # Initialize the lists for holding the logs, losses and accuracies
        self.losses = []
        self.acc = []
        self.val_losses = []
        self.val_acc = []
        self.logs = []
        self.x = []
        self.i = 1

    def on_epoch_end(self, epoch, logs=None):
        self.logs.append(logs)
        self.x.append(self.i)
        self.acc.append(logs.get('accuracy'))
        self.losses.append(logs.get('losses'))
        self.val_acc.append(logs.get('val_accuracy'))
        self.val_losses.append(logs.get('val_losses'))
        self.i += 1
        # If running in debug mode pop up a plot at each step
        if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
            plt.plot(self.x, self.acc, label="training", linestyle='-')
            plt.plot(self.x, self.val_acc, label="validation", linestyle='--')
            plt.legend()
            plt.show()


def load_training_and_validation_sets(image_dir="./data",
                                      planted_list_dir="./data",
                                      test_fraction=0.3,
                                      size=CUTOUT_DIMENSION,
                                      random=True,
                                      num_samples=1000,
                                      num_per_pair=2,
                                      num_pairs=None):
    """
    Loads from disk the image cutout sections (some with moving sources, some without) and returns as
    two groups of inputs shaped for the CNN model we will use.

    :type num_per_pair: int
    :param image_dir: directory containing image patches that will be loaded from disk
    :param planted_list_dir: directory containing files with lists of planted sources.
    :param test_fraction: what fraction of the input set will be used for training.
    :param size: what size of cutouts should we generate from the input images.
    :param random: should the grid of cutouts be regular (False) or randomized (True)
    :param num_samples: How many cutout samples should we gather.
    :param num_pairs: How many image pairs to use, None ==> All possible.
    :return: list of training and validation arrays.
    :rtype: np.array, np.array, nd.array, nd.array
    """

    # Load the image data from disk
    image_pairs = data_model.build_image_pair_list(image_directory=image_dir, num_pairs=num_pairs,
                                                   num_per_pair=num_per_pair)
    table_of_planted_sources = data_model.build_table_of_planted_sources(plant_list_directory=planted_list_dir)
    source_cutouts, source_targets, blank_cutouts = data_model.cut(image_pairs,
                                                                   table_of_planted_sources,
                                                                   size=size,
                                                                   random=random,
                                                                   num_samples=num_samples)

    # the data_model.cut has a pair index first index and we don't care about image pair
    # index here.  Reshape the data arrays to remove grouping by image pairs
    # Need to use the same number of blanks as source cutouts.
    logging.debug("Loading data resulted in {} possible source cutouts.".format(source_cutouts.shape[0]))
    logging.debug("Loading data resulted in {} possible blank cutouts.".format(blank_cutouts.shape[0]))
    num_of_sets_to_use = min(blank_cutouts.shape[0], source_cutouts.shape[0])
    logging.debug("This is how many cutouts we will use: {}".format(num_of_sets_to_use))
    if num_of_sets_to_use < 1:
        raise ValueError("Too little data for training.")
    # image_cutouts_grouped_by_pair hold the image cutouts in an array with each element having shape 2xndimx2ndim
    np.random.shuffle(source_cutouts)
    np.random.shuffle(blank_cutouts)
    image_cutouts = np.append(source_cutouts[0:num_of_sets_to_use],
                              blank_cutouts[0:num_of_sets_to_use], axis=0)
    # tar_bin holds 1 for the source and 0 for the blank cutouts.
    tar_bin = np.append(np.ones(num_of_sets_to_use), np.zeros(num_of_sets_to_use))
    return md.train_test_split(image_cutouts, tar_bin, test_size=test_fraction)


def get_cnn_model(channels=NUM_CHANNELS, dimension=CUTOUT_DIMENSION, kernel_size=2):
    """
    Get the CNN model used for the moving object detection process.
    :param channels: number of image cutouts packaged together, Currently we do image pairs.
    :param dimension: size, in pixels, of the image sections that will be provided.
    :param kernel_size: size, in pixels, of the Convolution kernel, 2 works well.
    # Not clear why kernel_size = 2 is used, what is the theory?
    :return: A keras CNN model.
    :rtype Model
    """

    # Tensorflow keep the channels AFTER the rows and columns.
    # We need to change it to channels first.
    keras.backend.set_image_data_format('channels_first')

    # input layer
    # The shape is (2,128,128). The first '2' means the number of channels.The second and the
    main_input1 = Input(shape=(channels, dimension, dimension), dtype='float32', name='main_input')
    logging.debug("Created Input object to build model around with shape: {}".format(main_input1.shape))
    # Conv2D is for the spatial convolution over images.
    # Usually better for images
    # MaxPooling2D is for calculating the maximum value for each patch of the feature map.
    # ReLu can be expressed as y = max(0,x). It returns the largest value, but ignore x below 0.
    # For classification, ReLu is mostly used hidden layers
    batch_normalize = False
    # Build a series of Conv2D layers starting with the main_input1 data.
    cnn = main_input1
    for node_count in [8, 16, 32, 64]:
        cnn = Conv2D(node_count, kernel_size=(kernel_size, kernel_size),
                     activation='relu', padding='same', strides=(1, 1))(cnn)
        # Pull the nodes
        cnn = MaxPooling2D(pool_size=(kernel_size, kernel_size))(cnn)

        # batch normalize on alternate layers.
        if batch_normalize:
            cnn = BatchNormalization()(cnn)
        else:
            batch_normalize = True

    # flatten the final Convolution layer to 1D.
    fc = Flatten(name='flat1')(cnn)

    # Dense is the regular densely connected NN layers
    # A dense layer represents a matrix vector multiplication.
    # The values in the matrix are the trainable parameters which get updated during backpropagation.
    # A dense layer thus is used to change the dimensions of your vector.
    # Mathematically speaking, it applies a rotation, scaling, translation transform to your vector.

    # A dropout layer is used for regularization where you randomly set some of the dimensions of
    # your input vector to be zero with probability keep_prob.
    # A dropout layer does not have any trainable parameters.
    # Dropout technique is to prevent over fitting.
    for node_count in [256, 64]:
        fc = Dense(node_count, activation='relu')(fc)
        fc = Dropout(0.25)(fc)

    # For the output layer, we can use either tanh or sigmoid.
    # The reason I use sigmoid (non-linear) function here is because it is a classification problem.
    # The sigmoid can be expressed as y = 1/(1+e^(-x))
    # We can set the rule for this activation function.
    # For example, if the y is less than 0.5, it outputs 0 and if the y is more than 0.5, it outputs 1.
    main_output = Dense(1, activation='sigmoid')(fc)
    model = Model(inputs=main_input1, outputs=main_output)

    # "compile" configures the model.
    # The first argument is an optimizer.
    # An optimizer is an algorithm to change weights and learning rate in order to reduce the losses.
    # The second argument is a loss. Machines learn by means of a loss function.
    # It evaluates how well specific algorithm models the given data.
    # The third argument is a metric. I used only accuracy here.
    # The definition of accuracy is given in the later cells.
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    # Display the summary of the model.
    logging.debug(model.summary())
    return model


def train_and_validate_the_model(model, training_data, training_classes, validation_data, validation_classes,
                                 batch_size=64, epochs=50):
    """
    This is the actual fitting part.  Train your dragon.

    :param model: the CNN model to be trained.
    :param training_data: list of image cutouts (grouped by image pairs) to train with.
    :param training_classes: list of category values (0 or 1) for training set, same order as training_data.
    :param validation_data: List of image cutouts (grouped by images pairs) to validate with.
    :param validation_classes: list of category values (0 or 1) for validation set, same order as validation_data
    :param batch_size: How many training sets to pass in a single batch (memory consideration)
    :param epochs: How many training loops to build model (time consideration and over fitting concern)
    :return: The history object used to examine how the training went.
    :rtype: Model.History
    """

    # validation_data = keras evaluates the loss and accuracy at the end of each epoch.
    # This data is not used for training the model.
    # shuffle = keras shuffles training inputs in batch-sized chunks.
    # callbacks = a function that is called at the end of each epoch.

    history = model.fit(training_data, training_classes,
                        validation_data=[validation_data, validation_classes],
                        shuffle=True,
                        epochs=epochs,
                        batch_size=batch_size,
                        callbacks=[TrainingPlot()])
    return history


def main(model_filename="trained_model.ker"):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num-samples', help='How many samples images to create per pair?', default=2000, type=int)
    parser.add_argument('--num-pairs', help='How many pairs to generate?', default=100, type=int)
    parser.add_argument('--num-per-pair', help='How many images in a pair?', default=NUM_CHANNELS, type=int)
    parser.add_argument('--random', help='Should we draw out random image sections or follow a grid?',
                        action='store_true')
    parser.add_argument('--epochs', help='How many training epochs to run?', default=50, type=int)
    parser.add_argument('plant_list_dir', help='Directory containing lists of planted sources in plantList format')
    parser.add_argument('image_dir', help='Directory containing images used to train the model.')
    parser.add_argument('--test-fraction', help='Fraction of the model used to validate with.', default=0.3,
                        type=float)
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    parser.add_argument('--cutout-dimension', help='size, in pixels, of cutouts dimension to work with', type=int,
                        default=CUTOUT_DIMENSION)
    args = parser.parse_args()
    log_levels = {'DEBUG': logging.DEBUG, 'INFO': logging.INFO, 'ERROR': logging.ERROR}
    logging.basicConfig(level=log_levels[args.log_level])
    logging.info("Building input dataset for training.")

    training_data, validation_data, training_classes, validation_classes = \
        load_training_and_validation_sets(num_samples=args.num_samples,
                                          num_pairs=args.num_pairs,
                                          image_dir=args.image_dir,
                                          random=args.random,
                                          test_fraction=args.test_fraction,
                                          planted_list_dir=args.plant_list_dir,
                                          num_per_pair=args.num_per_pair,
                                          size=args.cutout_dimension)
    logging.info("Constructing the model framework.")
    model = get_cnn_model(channels=args.num_per_pair, dimension=args.cutout_dimension)
    logging.info("Training the model.")
    history = train_and_validate_the_model(model, training_data, training_classes,
                                           validation_data, validation_classes, epochs=args.epochs)
    logging.info("Plotting the history.")
    plot_training_outcome(history.history, output_file_base='grid')
    model.save(f'{args.cutout_dimension}_{args.num_per_pair}_{model_filename}')
    return


if __name__ == '__main__':
    sys.exit(main())
