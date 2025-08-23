# Maintainer: Chan Beom Park <cbpark@gmail.com>

pkgname=yam2-git
_pkgname=YAM2
pkgdesc="Yet another library for the M2 variables"
url=https://github.com/cbpark/${_pkgname}
pkgver=2.0.0
pkgrel=1
arch=('x86_64')
license=('BSD')
makedepends=('git' 'cmake')
depends=('nlopt')
optdepends=('root')
source=("git+${url}.git")
md5sums=('SKIP')

build() {
    cd "${_pkgname}"/build

    cmake .. \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/usr

    make
}

package() {
    cd "${_pkgname}"/build

    make DESTDIR="${pkgdir}" install

    mkdir -p "${pkgdir}/usr/share/license/${pkgname}"
    cd ..
    install -Dm 644 LICENSE -t "${pkgdir}/usr/share/license/${pkgname}/LICENSE"
}
