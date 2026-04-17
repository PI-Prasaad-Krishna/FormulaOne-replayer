import { useEffect, useState } from 'react';

const releaseBaseUrl = 'https://github.com/PI-Prasaad-Krishna/FormulaOne-replayer/releases/latest';

const downloads = [
  {
    name: 'Windows',
    filename: 'F1 Visualizer Setup.msi',
    details: 'Recommended for Windows 10 and 11. Fast install, native desktop shortcuts, ready to run.',
    badge: 'Most downloaded',
    accent: 'linear-gradient(135deg, rgba(255, 59, 48, 0.24), rgba(255, 255, 255, 0.08))',
    href: releaseBaseUrl,
  },
  {
    name: 'macOS',
    filename: 'F1 Visualizer.dmg',
    details: 'Built for Mac users who want a clean install and smooth launch on Apple Silicon and Intel.',
    badge: 'macOS ready',
    accent: 'linear-gradient(135deg, rgba(255, 255, 255, 0.18), rgba(255, 255, 255, 0.05))',
    href: releaseBaseUrl,
  },
  {
    name: 'Linux',
    filename: 'F1 Visualizer.AppImage',
    details: 'A portable Linux build for quick setup across supported distributions.',
    badge: 'Portable',
    accent: 'linear-gradient(135deg, rgba(0, 194, 255, 0.22), rgba(255, 255, 255, 0.06))',
    href: releaseBaseUrl,
  },
];

const highlights = [
  {
    title: 'Latest release',
    value: 'Available now',
  },
  {
    title: 'Platforms supported',
    value: 'Windows, macOS, Linux',
  },
  {
    title: 'Install experience',
    value: 'Fast and polished',
  },
];

const simpleSteps = [
  {
    title: 'Choose your platform',
    text: 'Pick Windows, macOS, or Linux from the download cards below.',
  },
  {
    title: 'Open the release page',
    text: 'Each button takes you to the official GitHub release for the latest build.',
  },
  {
    title: 'Install when ready',
    text: 'The installer files will be published there as soon as they are available.',
  },
];

function App() {
  const [mounted, setMounted] = useState(false);

  useEffect(() => {
    setMounted(true);

    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add('is-visible');
          }
        });
      },
      {
        threshold: 0.2,
        rootMargin: '0px 0px -10% 0px',
      },
    );

    document.querySelectorAll('[data-reveal]').forEach((node) => observer.observe(node));

    return () => observer.disconnect();
  }, []);

  return (
    <div className={`site ${mounted ? 'is-ready' : ''}`}>
      <div className="ambient ambient-left" />
      <div className="ambient ambient-right" />
      <div className="grid-overlay" />

      <header className="topbar">
        <a className="brand" href="#top" aria-label="F1 Visualizer home">
          <span className="brand-mark" aria-hidden="true">
            <svg viewBox="0 0 48 48" role="img" aria-hidden="true">
              <path d="M7 27h18l-3 7h9l10-17h-7l4-8H21l-8 10H7z" />
            </svg>
          </span>
          <span>
            <strong>F1 Visualizer</strong>
            <small>Official downloads</small>
          </span>
        </a>

        <nav className="nav">
          <a href="#downloads">Download</a>
          <a href="#support">Support</a>
        </nav>

        <a className="release-pill" href={releaseBaseUrl} target="_blank" rel="noreferrer">
          Latest release
        </a>
      </header>

      <main id="top">
        <section className="hero">
          <div className="hero-copy" data-reveal>
            <div className="eyebrow">Official release download</div>
            <h1>Download F1 Visualizer for your desktop.</h1>
            <p>
              Get the latest version for Windows, macOS, or Linux. Clean installation, premium
              presentation, and a release built to be easy to trust.
            </p>

            <div className="cta-row">
              <a className="primary-cta" href="#downloads">
                Download latest release
              </a>
              <a className="secondary-cta" href={releaseBaseUrl} target="_blank" rel="noreferrer">
                View on GitHub Releases
              </a>
            </div>

            <div className="stat-row">
              {highlights.map((item) => (
                <div key={item.title} className="stat-card">
                  <span>{item.title}</span>
                  <strong>{item.value}</strong>
                </div>
              ))}
            </div>
          </div>

          <div className="hero-visual" data-reveal>
            <div className="track-panel vehicle-panel">
              <div className="track-label">Formula 1-inspired styling</div>
              <div className="hero-accent-line" aria-hidden="true" />
              <svg className="car-illustration" viewBox="0 0 720 420" role="img" aria-label="Formula One style race car illustration">
                <defs>
                  <linearGradient id="carBody" x1="0" x2="1" y1="0" y2="1">
                    <stop offset="0%" stopColor="#ffffff" />
                    <stop offset="42%" stopColor="#f4f4f4" />
                    <stop offset="100%" stopColor="#b8b8b8" />
                  </linearGradient>
                  <linearGradient id="carRed" x1="0" x2="1" y1="0" y2="1">
                    <stop offset="0%" stopColor="#ff5a4f" />
                    <stop offset="100%" stopColor="#b91208" />
                  </linearGradient>
                </defs>
                <g opacity="0.22">
                  <path d="M100 308c124-26 252-26 420 0" fill="none" stroke="#ffffff" strokeWidth="6" strokeLinecap="round" />
                  <path d="M88 232c153-68 363-68 532 0" fill="none" stroke="#ffffff" strokeOpacity="0.16" strokeWidth="8" strokeLinecap="round" />
                </g>
                <g className="car-motion car-motion-one">
                  <path d="M186 240c-10-16 7-33 32-33h104l17-23h74c26 0 52 13 65 32l10 15h-54l-11 16h-78l-15 19H250c-23 0-50-10-64-26Z" fill="url(#carBody)" />
                  <path d="M270 188h130l27 38h-84l-13-18h-61z" fill="url(#carRed)" />
                  <path d="M169 250h61l-18 26h-34c-16 0-26-10-26-24 0-1 0-2 1-2 3 0 9 0 16 0Z" fill="#ffffff" />
                  <path d="M456 238h64l15 18h-36l-17-20z" fill="#ffffff" />
                  <circle cx="230" cy="270" r="39" fill="#111318" />
                  <circle cx="230" cy="270" r="20" fill="#f7f7f7" />
                  <circle cx="475" cy="270" r="39" fill="#111318" />
                  <circle cx="475" cy="270" r="20" fill="#f7f7f7" />
                  <circle cx="336" cy="225" r="10" fill="#ffffff" opacity="0.9" />
                </g>
                <g className="car-motion car-motion-two">
                  <path d="M218 154h86l17-13h82c18 0 34 11 40 26l2 6h-44l-8 10h-60l-14 15h-71c-18 0-31-13-31-28 0-7 1-13 1-16Z" fill="url(#carBody)" opacity="0.82" />
                  <path d="M286 138h90l15 24h-64l-9-12h-40z" fill="url(#carRed)" opacity="0.85" />
                  <circle cx="248" cy="196" r="22" fill="#111318" opacity="0.85" />
                  <circle cx="420" cy="196" r="22" fill="#111318" opacity="0.85" />
                </g>
              </svg>

              <div className="vehicle-livery" aria-hidden="true">
                <span />
                <span />
                <span />
              </div>

              <div className="build-card">
                <div>
                  <span className="muted">Download destination</span>
                  <strong>GitHub Releases</strong>
                </div>
                <div className="build-meta">
                  <span>Latest version</span>
                  <span>Trusted source</span>
                </div>
              </div>
            </div>
          </div>
        </section>

        <section className="downloads-section" id="downloads">
          <div className="section-heading" data-reveal>
            <span>Downloads</span>
            <h2>Select the installer for your device.</h2>
          </div>

          <div className="download-grid">
            {downloads.map((item) => (
              <a
                key={item.name}
                className="download-card"
                href={item.href}
                target="_blank"
                rel="noreferrer"
                style={{ background: item.accent }}
                data-reveal
              >
                <div className="download-card-top">
                  <span className="download-badge">{item.badge}</span>
                  <span className="download-platform">{item.name}</span>
                </div>

                <div className="download-main">
                  <strong>{item.filename}</strong>
                  <p>{item.details}</p>
                </div>

                <div className="download-footer">
                  <span>Open latest release</span>
                  <span aria-hidden="true">
                    <svg viewBox="0 0 24 24">
                      <path d="M13 5l7 7-7 7M20 12H4" />
                    </svg>
                  </span>
                </div>
              </a>
            ))}
          </div>
        </section>

        <section className="support-panel" id="support" data-reveal>
          <div>
            <span className="eyebrow">Need help</span>
            <h2>Choose the platform you use and open the latest release.</h2>
          </div>

          <div className="step-grid">
            {simpleSteps.map((step) => (
              <div key={step.title} className="step-card">
                <strong>{step.title}</strong>
                <p>{step.text}</p>
              </div>
            ))}
          </div>
        </section>
      </main>
    </div>
  );
}

export default App;