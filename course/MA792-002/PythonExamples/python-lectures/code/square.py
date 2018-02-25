from image import Image

class SquareImage(Image):
    def __init__(self, side):
        super(SquareImage, self).__init__(side, side)

square = SquareImage()
print square.dimensions()
