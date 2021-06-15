from unittest import TestCase
from . import sns
import numpy


class Test(TestCase):

    def test_weighted_quantile(self):
        n = numpy.random.choice(numpy.arange(50), [50, ], replace=False)
        image_stack = numpy.array([numpy.ones((10, 10))*i for i in n])
        image_weights = 0*image_stack+1/50.
        wq = sns.weighted_quantile(image_stack, 0.50001, image_weights)
        self.assertEqual(wq.shape[0], image_stack.shape[1])
        self.assertEqual(wq.shape[1], image_stack.shape[2])
        self.assertAlmostEqual(wq[5, 5], 25, 2)
