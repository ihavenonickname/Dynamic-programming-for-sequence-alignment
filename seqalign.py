class SequencesAligner:
  def __init__(self, seq1, seq2, match=1, mismatch=-1, gap=-1, local=False):
    self.match = match
    self.mismatch = mismatch
    self.gap = gap
    self.local = local
    self.seq1 = seq1 # This guy represents columns in matrix.
    self.seq2 = seq2 # This one represents lines.

  def initmatrix(self):
    # Matrix is an array of arrays, each inner array is a line.
    self.matrix = []

    columncounter = xrange(len(self.seq1) + 1)

    # First line is a 0 followed by gap times i, being i the index of the i-th
    # column.
    firstline = [self.verifyvalue(self.gap * i) for i in columncounter]

    self.matrix.append(firstline)

    # After the first one, each line starts with a cell being gap times i.
    for i in xrange(1, len(self.seq2) + 1):
      self.matrix.append([self.verifyvalue(i * self.gap)])

  def neighbours(self, column, line):
    up = self.matrix[line - 1][column] # Above.
    left = self.matrix[line][column - 1] # At left.
    diagn = self.matrix[line - 1][column - 1] # At the left of the above one.

    return up, left, diagn

  def hasmatch(self, column, line):
    # Is the column-th item in seq1 equal to line-th item in seq2?
    return self.seq1[column - 1] == self.seq2[line - 1]

  def thegreater(self, a, b, c, index=False):
    # List of 2-ples in the form (value, index), for preserving the original
    # position on the list after sorting it.
    values = [(a, 1), (b, 2), (c, 3)]

    # Last value after sorting the list in asc order.
    greater = sorted(values, key=lambda x: x[0])[2]

    if index:
      return greater[1]

    return greater[0]

  def verifyvalue(self, value):
    # If the alignment is local, negative values are considered zero.
    if self.local and value < 0:
      return 0

    return value

  def calcneighbours(self, up, left, diagn, hasmatch):
    # 1. The value coming from above is the value of the cell above plus
    # the value of gap.
    #
    # 2. The value coming from left is the value of the cell at left plus
    # the value of gap.
    #
    # 3. The value coming from diagonal is the value of the cell in
    # upper-left diagonal plus the value of match or mismatch.

    up = up + self.gap
    left = left + self.gap

    if hasmatch:
      diagn += self.match
    else:
      diagn += self.mismatch

    return up, left, diagn

  def generatevalue(self, column, line):
    # To generate the value of a given cell, its neighbours have to be
    # already generated, because the value of the actual cell dependes of
    # the value of its neighbours after applying the algorithm described
    # in calcneighbours function.
    # The greater value of these values is the value given to the actual cell.

    up, left, diagn = self.neighbours(column, line)

    hasmatch = self.hasmatch(column, line)

    up, left, diagn = self.calcneighbours(up, left, diagn, hasmatch)

    return self.thegreater(up, left, diagn)

  def mountmatrix(self):
    columncounter = xrange(1, len(self.seq1) + 1)
    linecounter = xrange(1, len(self.seq2) + 1)

    # Generates and appends missing values of each cell. The first line is not
    # included because it was filled in initmatrix. Same for the first column.
    for lineindex in linecounter:
      for columnindex in columncounter:
        value = self.generatevalue(columnindex, lineindex)
        self.matrix[lineindex].append(self.verifyvalue(value))

    return self.matrix

  def gloriousroute(self):
    # In this aproach, the route is got by finding which cell generated
    # another, starting from the last cell generated (the one at bottom-right
    # corner in matrix).
    route = []

    # Bottom-right cell.
    column = len(self.seq1)
    line = len(self.seq2)

    # While scanning valid indexes.
    while column >= 0 and line >= 0:
      # Add current cell to the list (a 2-ple).
      route.append((column, line))

      up, left, diagn = self.neighbours(column, line)

      hasmatch = self.hasmatch(column, line)

      up, left, diagn = self.calcneighbours(up, left, diagn, hasmatch)

      # Find out which cell generated current cell.
      nextcell = self.thegreater(up, left, diagn, index=True)

      # Switch to this cell.
      if nextcell == 1: # Up
        line -= 1
      elif nextcell == 2: # Left
        column -= 1
      else: # Diagn
        line -= 1
        column -= 1

    return route

  def alignment(self):
    newseq1 = ""
    newseq2 = ""

    # Route is traveled from back to front because new sequences are mounted
    # this way.
    route = self.gloriousroute()[::-1]

    i = 0
    seq1counter = 0
    seq2counter = 0

    while i + 1 < len(route):
      # If the column of current cell is different from the column of next
      # cell, then there was a move to left.
      fromleft = route[i][0] != route[i + 1][0]

      # If the line of current cell is different from the line of next cell,
      # then there was a move from above.
      fromabove = route[i][1] != route[i + 1][1]

      # If current cell came from left and from above, then it came from
      # upper-left diagonal, so both sequences are equal and have to be present
      # in the alignment.
      if fromabove and fromleft:
        newseq1 += self.seq1[seq1counter]
        newseq2 += self.seq2[seq2counter]

        seq1counter += 1
        seq2counter += 1

      # If current cell came only from above, then there is an item in the
      # second sequence that wasn't present in the first one, so an underline
      # is used to present it.
      elif fromabove:
        newseq1 += "_"
        newseq2 += self.seq2[seq2counter]

        seq2counter += 1

      # If current cell came only from left, then there is an item in the first
      # sequence that wasn't present in the second one, so an underline is used
      # to present it.
      else:
        newseq1 += self.seq1[seq1counter]
        newseq2 += "_"

        seq1counter += 1

      i += 1

    return newseq1, newseq2