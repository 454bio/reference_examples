import cv2

def findCircles2(img, minr=40, maxr=70):
    grayhighlow = img.copy()
    grayhighlow[grayhighlow < 22] = 0
    grayhighlow[grayhighlow > 0] = 255

    w2 = int(grayhighlow.shape[1] / 10)
    h2 = int(grayhighlow.shape[0] / 10)
    img2 = cv2.resize(grayhighlow, (w2,h2), interpolation = cv2.INTER_AREA)
    spots = []
    blur = cv2.GaussianBlur(img2,(15,15),cv2.BORDER_DEFAULT)
    h = grayhighlow.shape[0]
    w = grayhighlow.shape[1]
    for y in range(1,h2-1):
        for x in range(1,w2-1):
            if blur[y,x] > blur[y-1,x] and blur[y,x] >= blur[y+1,x] and blur[y,x] > blur[y,x-1] and blur[y,x] >= blur[y,x+1]:
                spots.append([x,y])
    print('found %d local max' % len(spots))
    '''
    if show_plots:
        fig, ax = plt.subplots()
        ax.imshow(blur, cmap='gray')
        for spot in spots:
            circle = plt.Circle((spot[0], spot[1]), 2, color='r')
            ax.add_patch(circle)
    '''
    spots2 = []
    for spot in spots:
        spot[0] = spot[0]*10 + 5
        spot[1] = spot[1]*10 + 5
        if spot[0] >= maxr and spot[0] < (w-maxr) and spot[1] >= maxr and spot[1] < (h-maxr):
            spots2.append([spot[0], spot[1]])
    spots = spots2

    newspots = []
    for spot in spots:
        bestr = 0
        for r in range(minr, maxr+1):
            xsum = 0
            ysum = 0
            for i in range(-r,r+1):
                xsum += grayhighlow[spot[1],spot[0]+i]
                ysum += grayhighlow[spot[1]+i,spot[0]]
            rsum = 0.95 * 255 * (2*r+1)
            if xsum >= rsum and ysum >= rsum:
                bestr = r
        if bestr > 0:
            newspots.append([spot[0], spot[1], bestr])
    print('found %d valid circles' % len(newspots))
    return newspots
